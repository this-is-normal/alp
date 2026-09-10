//! **Quasar**: a quasi-abelian code commitment switched into BaseFold — the
//! fast-prover polylog-verifier backend.
//!
//! Both the QA code and the foldable code are ChaCha-seeded, so an x86 prover
//! and a riscv32 verifier build identical public parameters.

use p3_basefold::{BasefoldPcs, RandomFoldableCode};
use p3_qa_code::QaCode;
use p3_quasar::{QuasarMlPcs, QuasarParams, QuasarPcs};
use p3_util::log2_ceil_usize;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::hash::{F, L, MyMmcs, build_mmcs};

/// Human-readable backend name, for measurement output.
pub const NAME: &str = "Quasar";

pub type CodePcs = QuasarMlPcs<F, MyMmcs, L>;

/// QA index `c` = inverse rate of the committing code.
const INDEX: usize = 2;
/// ed25519 is 255 bits — the field size the QA distance bound is taken at.
const FIELD_BITS: usize = 255;
/// Targets for the column-sampling error and for the distance bound's own
/// failure probability (the paper's Table 4 uses 103 for both).
const SECURITY_BITS: usize = 103;
const DISTANCE_FAILURE_BITS: usize = 103;
/// Foldable-code inverse rate for the code-switching openings, and its
/// list-decoding proximity radius (1/4 ⇒ δ ≈ 0.215, 1/16 ⇒ δ ≈ 0.350).
///
/// **MEASURED: keep 1/4 here, unlike the standalone BaseFold backend.** 1/16
/// buys 241 → 136 queries, but Quasar's recursion guest rebuilds this code
/// eagerly, and a 4× longer codeword means 4× the twist sampling: the OpenVM
/// guest went 1.26B → 2.13B instructions at 16 signatures. The rate is only
/// worth lowering here once the inner code is built lazily.
const BF_INV_RATE: usize = 4;
const BF_RADIUS: f64 = 0.215;
/// Proof-of-work bits ground before **each** query phase — the committed
/// columns, and the inner Basefold queries.
///
/// Each bit buys a bit of query soundness, so both counts are derived from the
/// reduced target `λ − POW_BITS` while the total stays the same. The prover pays
/// ~2^POW_BITS hashes per grind, and there are several grinds per proof (the
/// column phase and each inner Basefold opening, for both committed traces).
/// 16 bits costs ~65k hashes per grind — milliseconds under the BabyBear hash —
/// for a ~16% cut in both query counts. 20 bits was measurably too slow under
/// the 255-bit ed25519 hash, where a permutation is ~100x dearer.
const POW_BITS: usize = 16;
/// ChaCha seeds for the two public codes.
const QA_CODE_SEED: u64 = 0x0A_C0DE;
const BF_CODE_SEED: u64 = 0xBA5E_F01D;

/// The resolved public shape for a trace, so a recursion host can ship it
/// instead of making the guest rediscover it.
///
/// [`QuasarParams::choose`] bisects on `f64` transcendentals
/// (`log`/`pow`/`expm1`). On a 64-bit host that is microseconds; on rv32im
/// there is no FPU, so every one of them is a softfloat call, and the search
/// measured **408 M instructions — 63 % of the whole OpenVM recursion guest**,
/// dwarfing the verification it was configuring. The shape is a public
/// parameter of a fixed trace geometry, exactly like the coefficient
/// commitment, so it belongs in the verifier's inputs, not in its work.
#[must_use]
pub fn resolved_params(log_height: usize, width: usize) -> (QuasarParams, usize) {
    (choose_params(log_height, width), basefold_queries())
}

/// Build one Quasar PCS sized for a `2^log_height`-row, `width`-column trace.
///
/// The row/column split comes from [`QuasarParams::choose`], which minimises
/// predicted proof size; the column count comes from the proven QA distance,
/// clamped to the codeword so tiny instances stay well-formed (and, as
/// `achieved_security_bits` reports, honestly under-secured).
fn choose_params(log_height: usize, width: usize) -> QuasarParams {
    let stacked = log_height + log2_ceil_usize(width);
    let column_target = SECURITY_BITS.saturating_sub(POW_BITS).max(1);
    let bf_queries = basefold_queries();
    let params = QuasarParams::choose(
        stacked,
        INDEX,
        FIELD_BITS,
        column_target,
        DISTANCE_FAILURE_BITS,
        BF_INV_RATE,
        bf_queries,
    )
    .or_else(|| {
        // No shape reaches the target (a small trace): take the smallest-proof
        // shape at whatever soundness the codeword can carry, and say so.
        QuasarParams::choose_best_effort(
            stacked,
            INDEX,
            FIELD_BITS,
            column_target,
            DISTANCE_FAILURE_BITS,
            BF_INV_RATE,
            bf_queries,
        )
    })
    .expect("ed25519 admits a QA distance bound");
    QuasarParams::new(
        params.log_rows,
        params.log_cols,
        INDEX,
        params.num_queries.min(params.codeword_len()),
    )
}

/// Build one Quasar PCS sized for a `2^log_height`-row, `width`-column trace.
pub fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    let params = choose_params(log_height, width);
    let bf_queries = basefold_queries();
    let log_cols = params.log_cols;
    let code = QaCode::<F>::new_random(
        log_cols,
        INDEX,
        &mut ChaCha8Rng::seed_from_u64(QA_CODE_SEED),
    );
    let basefold = BasefoldPcs::new_with_grinding(
        RandomFoldableCode::<F>::new(log_cols, BF_INV_RATE, BF_CODE_SEED),
        build_mmcs(),
        bf_queries.min(BF_INV_RATE << log_cols),
        POW_BITS,
    );

    tracing::info!(
        "quasar c={} N1=2^{} N2=2^{}: {} of {} columns queried + {POW_BITS} pow bits \
         = {:.1}-bit sampling soundness; basefold {} queries at rate 1/{}",
        INDEX,
        params.log_rows,
        log_cols,
        params.num_queries,
        params.codeword_len(),
        params
            .achieved_security_bits(FIELD_BITS, DISTANCE_FAILURE_BITS)
            .unwrap_or(0.0)
            + POW_BITS as f64,
        bf_queries.min(BF_INV_RATE << log_cols),
        BF_INV_RATE,
    );

    let (pcs, key) = QuasarPcs::new(params, code, build_mmcs(), basefold);
    CodePcs::new(pcs.with_grinding(POW_BITS), key)
}

/// Basefold query repetitions for a 100-bit target, net of the grinding bits:
/// `q = ⌈(100 − POW_BITS) / −log₂(1 − δ)⌉`.
fn basefold_queries() -> usize {
    let per_query = -(1.0 - BF_RADIUS).log2();
    ((100.0 - POW_BITS as f64).max(1.0) / per_query).ceil() as usize
}
