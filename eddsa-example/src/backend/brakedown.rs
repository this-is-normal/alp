//! The recursive **Brakedown** linear-time code in a tensor PCS: linear-time
//! encoding, `O(√N)` proof and verifier.

use p3_brakedown::{BrakedownCode, BrakedownMlPcs, LinearCode, TensorPcs};
use p3_util::log2_ceil_usize;
#[cfg(feature = "bb-hash")]
use rand_chacha::ChaCha8Rng;
#[cfg(not(feature = "bb-hash"))]
use rand::rngs::SmallRng;
use rand::SeedableRng;

use crate::hash::{EF, F, L, MyMmcs, build_mmcs};

/// The tensor PCS's linear code.
pub type Code = BrakedownCode<F, EF>;
/// The multilinear PCS used by multi-stark.
/// Human-readable backend name, for measurement output.
pub const NAME: &str = "Brakedown";

pub type CodePcs = BrakedownMlPcs<F, EF, Code, MyMmcs, L>;

/// Fixed seed for the code's random expander matrices, so prover and verifier
/// build the identical linear code.
const CODE_SEED: u64 = 0xB2A4E_D0;
/// Row/column split strategy for the tensor commitment.
///
/// `MinSize` reproduces the old `optimal_log_rows` behaviour (smallest proof).
/// `MinArith`/`MinHashes`/`MinVerify` trade proof size for fewer verifier
/// operations — the axis that matters when a zkVM guest does the verifying.
const SPLIT: p3_brakedown::params::SplitStrategy = p3_brakedown::params::SplitStrategy::MinSize;
/// `false` = provable δ/3 radius (Ligero/Roth–Zémor); `true` = conjectured δ/2
/// (Diamond–Posen Conjecture-1).
const CONJECTURED: bool = false;

/// Build one tensor PCS sized for a `2^log_height`-row, `width`-column trace,
/// targeting `security_bits` of column-query soundness.
///
/// The committed polynomial stacks every trace column, so the scheme is sized
/// for the stacked variable count `log_height + ceil(log2(width))`, split by
/// [`SPLIT`].
pub fn build_pcs(log_height: usize, width: usize, security_bits: usize) -> CodePcs {
    let stacked = log_height + log2_ceil_usize(width);

    // Small-proof preset (Haböck SdigCode6, δ=0.070, weights 9/12) with a
    // distance-derived query count.
    let preset = p3_brakedown::params::small_proof();
    let t = preset.num_queries(security_bits, CONJECTURED);
    // The cost model the split is optimised against: BabyBear hashing when the
    // guest verifies, Poseidon2-over-ed25519 otherwise.
    #[cfg(feature = "bb-hash")]
    let model = p3_brakedown::params::SplitModel::bb_hash();
    #[cfg(not(feature = "bb-hash"))]
    let model = p3_brakedown::params::SplitModel::native();
    let log_rows = p3_brakedown::params::choose_log_rows(stacked, t, &preset, &model, SPLIT);
    let cols = 1 << (stacked - log_rows);

    // ChaCha under `bb-hash` (portable, so an x86 prover and a riscv32 verifier
    // build identical matrices); SmallRng otherwise, matching the older runs.
    #[cfg(feature = "bb-hash")]
    let mut code_rng = ChaCha8Rng::seed_from_u64(CODE_SEED);
    #[cfg(not(feature = "bb-hash"))]
    let mut code_rng = SmallRng::seed_from_u64(CODE_SEED);
    let code = Code::from_message_len_with_params(cols, &preset.params, &mut code_rng);

    let num_queries = core::cmp::min(t, code.codeword_len());
    CodePcs::new(TensorPcs::new(code, build_mmcs(), log_rows, num_queries))
}
