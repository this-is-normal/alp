//! Multi-stark backend for the EdDSA AIR, over the ed25519 base field.
//!
//! `EdDsaAir` declares a **preprocessed trace** (the fixed-base doubling ladder +
//! row selectors), so the config provides a *second* commitment scheme via
//! [`MultiStarkConfig::preprocessed_pcs`]; the preprocessed trace stacks a
//! different column count than the main trace, giving it a different stacked
//! arity, so it needs its own PCS instance sized for that arity.
//!
//! **Hash field (feature `bb-hash`).** By default the Merkle commitment and
//! Fiat-Shamir transcript hash with a Poseidon2 permutation over **ed25519**
//! (255-bit). With `--features bb-hash` the hash instead runs over **BabyBear**
//! (31-bit): the committed data, linear code, and sumcheck stay over ed25519, and
//! only the hash field changes — each ed25519 value is decomposed into
//! base-`|BabyBear|` limbs for hashing (see [`SplitFieldHasher`]) and the
//! challenger sponges over BabyBear while presenting an ed25519 interface (see
//! [`DuplexChallengerCF`]). This makes every Poseidon2 permutation ~100x cheaper
//! inside the 32-bit OpenVM verifier; the default ed25519 hash is kept for
//! host-side comparison and op-count calibration.

#[cfg(all(not(feature = "ligero"), not(feature = "basefold")))]
use p3_brakedown::BrakedownCode as RealBrakedownCode;
#[cfg(all(feature = "ligero", not(feature = "basefold")))]
use p3_brakedown::ReedSolomonCode;
#[cfg(not(feature = "basefold"))]
use p3_brakedown::{BrakedownMlPcs, LinearCode, TensorPcs};
#[cfg(feature = "basefold")]
use p3_basefold::{BasefoldMlPcs, BasefoldPcs, RandomFoldableCode};
use p3_ed25519::Ed25519;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multi_stark::config::MultiStarkConfig;
use p3_sumcheck::layout::{Layout, PrefixProver, Table, Witness};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_util::log2_ceil_usize;
// Needed for `seed_from_u64` on the ed25519 perm RNG (default) and the Brakedown code RNG
// (non-basefold); the `bb-hash`+`basefold` combo seeds neither in this module.
#[cfg(not(all(feature = "bb-hash", feature = "basefold")))]
use rand::SeedableRng;

// ---- ed25519-hash-only imports (default) ----
#[cfg(not(feature = "bb-hash"))]
use p3_challenger::DuplexChallengerPF;
#[cfg(not(feature = "bb-hash"))]
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Poseidon2Ed25519,
};
#[cfg(not(feature = "bb-hash"))]
use rand::rngs::SmallRng;

// ---- bb-hash-only imports ----
#[cfg(feature = "bb-hash")]
use p3_baby_bear::{BabyBear, Poseidon2BabyBear, default_babybear_poseidon2_16};
#[cfg(feature = "bb-hash")]
use p3_challenger::DuplexChallengerCF;
#[cfg(feature = "bb-hash")]
use p3_symmetric::SplitFieldHasher;
// Only the non-basefold Brakedown code path seeds a portable ChaCha RNG.
#[cfg(all(feature = "bb-hash", not(feature = "basefold")))]
use rand_chacha::ChaCha8Rng;

// ---- count-ops instrumentation (either hash field) ----
#[cfg(feature = "count-ops")]
use p3_symmetric::{CryptographicPermutation, Permutation};

/// Base (trace) field.
pub type F = Ed25519;
/// Challenge field: the degree-1 self-extension of ed25519 (`|EF| ≈ 2^255`).
pub type EF = Ed25519;

// =====================================================================================
// Hash stack A: Poseidon2-over-ed25519 (default).
// =====================================================================================

/// Poseidon2 state width (ed25519 Poseidon2 supports width 3).
#[cfg(not(feature = "bb-hash"))]
const WIDTH: usize = 3;
/// Sponge rate (challenger + hash); capacity = 1 (~255 bits).
#[cfg(not(feature = "bb-hash"))]
const RATE: usize = 2;
/// Merkle digest / sponge output size in field elements.
#[cfg(not(feature = "bb-hash"))]
const OUT: usize = 1;

#[cfg(all(not(feature = "bb-hash"), not(feature = "count-ops")))]
type Perm = Poseidon2Ed25519<WIDTH>;
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
type Perm = CountingPerm;

/// Counts every Poseidon2-over-Ed25519 permutation (a hash invocation). Locally the hash is over
/// Ed25519, so its muls are also counted in `MUL`; counting perms lets us subtract them to isolate
/// the PCS *arithmetic* muls (with `bb-hash` the hash is BabyBear, so this subtraction is moot).
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
#[derive(Clone, Debug)]
pub struct CountingPerm(Poseidon2Ed25519<WIDTH>);
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
impl Permutation<[Ed25519; WIDTH]> for CountingPerm {
    fn permute_mut(&self, input: &mut [Ed25519; WIDTH]) {
        p3_ed25519::counters::on_perm();
        self.0.permute_mut(input);
    }
}
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
impl CryptographicPermutation<[Ed25519; WIDTH]> for CountingPerm {}

#[cfg(not(feature = "bb-hash"))]
type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
#[cfg(not(feature = "bb-hash"))]
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
#[cfg(not(feature = "bb-hash"))]
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;

/// The general-`PrimeField` duplex challenger (ed25519 is not `PrimeField64`).
#[cfg(not(feature = "bb-hash"))]
pub type MyChallenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;

/// Fixed RNG seed for the (DEMO-grade) Poseidon2 round constants; seeding
/// deterministically makes the prover and verifier build the identical permutation.
#[cfg(not(feature = "bb-hash"))]
const POSEIDON2_SEED: u64 = 0xED25519_C0FFEE;

/// Build the fixed Poseidon2 permutation (ed25519 is not `PrimeField64`, so the
/// generic `new_from_rng` is used).
#[cfg(not(feature = "bb-hash"))]
fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(POSEIDON2_SEED);
    let inner = Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    );
    #[cfg(not(feature = "count-ops"))]
    {
        inner
    }
    #[cfg(feature = "count-ops")]
    {
        CountingPerm(inner)
    }
}

/// The Merkle MMCS (leaf hash + node compression) over the ed25519 permutation.
#[cfg(not(feature = "bb-hash"))]
fn build_mmcs() -> MyMmcs {
    MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0)
}

/// Calibration: ed25519 muls consumed by ONE Poseidon2 permutation (to subtract hashing from the
/// verify mul-count and isolate PCS arithmetic). Only meaningful for the ed25519 hash.
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
pub fn muls_per_perm() -> u64 {
    use p3_field::PrimeCharacteristicRing;
    let p = perm();
    let mut st = [F::ZERO; WIDTH];
    p3_ed25519::counters::reset();
    p.permute_mut(&mut st);
    p3_ed25519::counters::snapshot().mul
}

/// A bare seeded challenger — identical on the prove and verify sides (no domain separator).
#[cfg(not(feature = "bb-hash"))]
pub fn challenger() -> MyChallenger {
    MyChallenger::new(perm())
}

// =====================================================================================
// Hash stack B: Poseidon2-over-BabyBear (feature `bb-hash`).
// =====================================================================================

/// BabyBear Poseidon2 state width / sponge rate / digest size (31-bit hash field).
#[cfg(feature = "bb-hash")]
const BB_WIDTH: usize = 16;
#[cfg(feature = "bb-hash")]
const BB_RATE: usize = 8;
/// Merkle digest = 8 BabyBear ≈ 248-bit (~124-bit collision resistance).
#[cfg(feature = "bb-hash")]
const BB_DIGEST: usize = 8;

#[cfg(all(feature = "bb-hash", not(feature = "count-ops")))]
type PermBB = Poseidon2BabyBear<BB_WIDTH>;
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
type PermBB = CountingPerm;

/// Wraps the BabyBear Poseidon2 permutation to count every invocation (hash permutation)
/// during profiling. Value-identical to the inner permutation.
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
#[derive(Clone, Debug)]
pub struct CountingPerm(Poseidon2BabyBear<BB_WIDTH>);
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
impl Permutation<[BabyBear; BB_WIDTH]> for CountingPerm {
    fn permute_mut(&self, input: &mut [BabyBear; BB_WIDTH]) {
        p3_ed25519::counters::on_hash();
        self.0.permute_mut(input);
    }
}
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
impl CryptographicPermutation<[BabyBear; BB_WIDTH]> for CountingPerm {}

/// Inner BabyBear→BabyBear sponge (hashes the limb stream).
#[cfg(feature = "bb-hash")]
type InnerHash = PaddingFreeSponge<PermBB, BB_WIDTH, BB_RATE, BB_DIGEST>;
/// Leaf hasher: ed25519 elements → BabyBear digest via base-|BabyBear| limbs.
#[cfg(feature = "bb-hash")]
type MyHash = SplitFieldHasher<Ed25519, BabyBear, InnerHash>;
#[cfg(feature = "bb-hash")]
type MyCompress = TruncatedPermutation<PermBB, 2, BB_DIGEST, BB_WIDTH>;
#[cfg(feature = "bb-hash")]
type MyMmcs = MerkleTreeMmcs<Ed25519, BabyBear, MyHash, MyCompress, 2, BB_DIGEST>;

/// Transcript challenger: sponges over BabyBear, observes/samples ed25519.
#[cfg(feature = "bb-hash")]
pub type MyChallenger = DuplexChallengerCF<Ed25519, BabyBear, PermBB, BB_WIDTH, BB_RATE>;

/// Build the fixed BabyBear Poseidon2 permutation.
#[cfg(feature = "bb-hash")]
fn perm_bb() -> PermBB {
    // OpenVM-matching Horizen RC16 constants (identical to the poseidon2 precompile).
    #[cfg(not(feature = "count-ops"))]
    {
        default_babybear_poseidon2_16()
    }
    #[cfg(feature = "count-ops")]
    {
        CountingPerm(default_babybear_poseidon2_16())
    }
}

/// The Merkle MMCS: SplitFieldHasher leaves + BabyBear-permutation node compression.
#[cfg(feature = "bb-hash")]
fn build_mmcs() -> MyMmcs {
    MyMmcs::new(MyHash::new(InnerHash::new(perm_bb())), MyCompress::new(perm_bb()), 0)
}

/// A bare seeded challenger — identical on the prove and verify sides (no domain separator).
#[cfg(feature = "bb-hash")]
pub fn challenger() -> MyChallenger {
    MyChallenger::new(perm_bb())
}

// =====================================================================================
// Shared: layout, code, PCS selection, parameters, and the multi-stark config.
// =====================================================================================

type L = PrefixProver<F, EF>;

/// The tensor PCS's linear code. Default = the recursive **Brakedown** code
/// (linear-time). With `--features ligero` it becomes a **Reed–Solomon** code
/// (high distance → few queries, but superlinear Vandermonde encoding).
#[cfg(all(not(feature = "ligero"), not(feature = "basefold")))]
pub type Code = RealBrakedownCode<F, EF>;
#[cfg(all(feature = "ligero", not(feature = "basefold")))]
pub type Code = ReedSolomonCode<F>;

/// The multilinear PCS used by multi-stark. Default = the Brakedown tensor PCS; with
/// `--features basefold` = the ed25519-native Basefold PCS (polylog verifier).
#[cfg(not(feature = "basefold"))]
pub type CodePcs = BrakedownMlPcs<F, EF, Code, MyMmcs, L>;
#[cfg(feature = "basefold")]
pub type CodePcs = BasefoldMlPcs<F, MyMmcs, L>;

/// Fixed seed for the Brakedown code's random expander matrices, so prover and
/// verifier build the identical linear code.
#[cfg(all(not(feature = "ligero"), not(feature = "basefold")))]
const BRAKEDOWN_CODE_SEED: u64 = 0xB2A4E_D0;
/// Reed–Solomon inverse rate `n/k` (distance δ = 1 − 1/rate). Rate 2 → δ = 0.5.
#[cfg(feature = "ligero")]
const LIGERO_RATE: usize = 2;
/// Default target soundness (bits) for `config_for_eddsa` — Brakedown query count and the Ligero
/// query count both derive from it. Override per-call via `config_for_eddsa_with_security`.
const TARGET_SECURITY_BITS: usize = 100;
/// `false` = provable δ/3 radius (Ligero/Roth–Zémor); `true` = conjectured δ/2
/// (Diamond–Posen Conjecture-1, the Reed–Solomon-proven analogue).
#[cfg(all(not(feature = "ligero"), not(feature = "basefold")))]
const BRAKEDOWN_CONJECTURED: bool = false;

/// Row/column split strategy for the tensor commitment (see `p3_brakedown::params`).
/// `MinSize` = smallest proof (the previous `optimal_log_rows` behaviour); `MinArith`/`MinHashes`/
/// `MinVerify` trade proof size for fewer verifier ops/hashes — relevant when the guest verifies.
#[cfg(all(not(feature = "ligero"), not(feature = "basefold")))]
const BRAKEDOWN_SPLIT: p3_brakedown::params::SplitStrategy = p3_brakedown::params::SplitStrategy::MinSize;

/// Number of column queries for a Reed–Solomon (Ligero) tensor test at
/// `security_bits` bits, under the **list-decoding-up-to-capacity** conjecture:
/// `t = ⌈ λ / log₂(rate) ⌉`. This is the aggressive (conjectural) choice
/// production STARKs use for concrete parameters; it is NOT unconditionally proven.
#[cfg(all(feature = "ligero", not(feature = "basefold")))]
fn ligero_num_queries(rate: usize, security_bits: usize) -> usize {
    let delta = 1.0 - 1.0 / rate as f64; // = 1 − ρ, list-decoding capacity radius
    let bits_per_query = -(1.0 - delta).log2(); // = log2(rate)
    (security_bits as f64 / bits_per_query).ceil() as usize
}
/// Padding floor. The tensor PCS commits at natural arity, so any positive arity is fine.
const BRAKEDOWN_FOLDING: usize = 1;

/// Build one Basefold PCS sized for a `2^log_height`-row, `width`-column trace.
/// (`security_bits` is Brakedown's knob; Basefold's aggressive q=288 is fixed here.)
#[cfg(feature = "basefold")]
fn build_pcs(log_height: usize, width: usize, _security_bits: usize) -> CodePcs {
    // Basefold commits the stacked polynomial (num_variables = log_height + ceil(log2 width)).
    let num_variables = log_height + log2_ceil_usize(width);
    let mmcs = build_mmcs();
    // Rate 1/4 foldable code, ChaCha-seeded (cross-platform).
    let code = RandomFoldableCode::<F>::new(num_variables, 4, 0xBA5E_F01D);
    // 100-bit AGGRESSIVE security (rate 1/4, list-decoding radius δ≈0.215):
    // q = ⌈100 / −log₂(1−δ)⌉ = 288, clamped to the codeword length (= inv_rate·2^num_variables)
    // so tiny instances stay well-formed.
    let num_queries = core::cmp::min(288, 4usize << num_variables);
    BasefoldMlPcs::new(BasefoldPcs::new(code, mmcs, num_queries))
}

/// Build one Brakedown/Ligero tensor PCS sized for a `2^log_height`-row, `width`-column trace.
///
/// The committed polynomial stacks every trace column, so the tensor scheme is
/// sized for the stacked variable count `log_height + ceil(log2(width))`. Brakedown picks its
/// row/column split via `choose_log_rows(BRAKEDOWN_SPLIT)`; the RS/Ligero path uses a square split.
#[cfg(not(feature = "basefold"))]
fn build_pcs(log_height: usize, width: usize, security_bits: usize) -> CodePcs {
    let stacked_num_variables = log_height + log2_ceil_usize(width);

    let mmcs = build_mmcs();

    // Small-proof Brakedown preset (Haböck SdigCode6, δ=0.070, weights 9/12) with a real
    // distance-derived query count `bd_t` (instead of the old demo constant).
    #[cfg(not(feature = "ligero"))]
    let bd_preset = p3_brakedown::params::small_proof();
    #[cfg(not(feature = "ligero"))]
    let bd_t = bd_preset.num_queries(security_bits, BRAKEDOWN_CONJECTURED);

    // Strategy-chosen row/column split for Brakedown; square for the RS/Ligero path.
    #[cfg(all(not(feature = "ligero"), feature = "bb-hash"))]
    let bd_model = p3_brakedown::params::SplitModel::bb_hash();
    #[cfg(all(not(feature = "ligero"), not(feature = "bb-hash")))]
    let bd_model = p3_brakedown::params::SplitModel::native();
    #[cfg(not(feature = "ligero"))]
    let log_rows = p3_brakedown::params::choose_log_rows(
        stacked_num_variables,
        bd_t,
        &bd_preset,
        &bd_model,
        BRAKEDOWN_SPLIT,
    );
    #[cfg(feature = "ligero")]
    let log_rows = stacked_num_variables / 2;
    let cols = 1 << (stacked_num_variables - log_rows);

    #[cfg(all(not(feature = "ligero"), not(feature = "bb-hash")))]
    let code = {
        // Recursive Brakedown code; SmallRng-seeded so prover and verifier build the identical code.
        let mut code_rng = SmallRng::seed_from_u64(BRAKEDOWN_CODE_SEED);
        RealBrakedownCode::from_message_len_with_params(cols, &bd_preset.params, &mut code_rng)
    };
    #[cfg(all(not(feature = "ligero"), feature = "bb-hash"))]
    let code = {
        // ChaCha (portable) so the x86 prover and riscv32 verifier build identical matrices.
        let mut code_rng = ChaCha8Rng::seed_from_u64(BRAKEDOWN_CODE_SEED);
        RealBrakedownCode::from_message_len_with_params(cols, &bd_preset.params, &mut code_rng)
    };
    #[cfg(feature = "ligero")]
    // Reed–Solomon `[n, k]` with k = cols (power of two) and n = rate·cols.
    let code = ReedSolomonCode::new(cols, LIGERO_RATE * cols);

    #[cfg(not(feature = "ligero"))]
    // Distance-derived (provable δ/3 by default). The chosen split MUST supply enough columns for
    // the full query count — never silently clamp to fewer queries, which would prove insecurely.
    let num_queries = {
        assert!(
            code.codeword_len() >= bd_t,
            "Brakedown split (cols={cols}) has codeword {} < the {bd_t} queries needed for {security_bits}-bit \
             security — insecure; increase the trace size or pick a higher-δ preset",
            code.codeword_len(),
        );
        bd_t
    };
    #[cfg(feature = "ligero")]
    // Security-derived: t = ⌈λ / −log₂(1 − δ)⌉, clamped to the codeword length.
    let num_queries = core::cmp::min(
        ligero_num_queries(LIGERO_RATE, TARGET_SECURITY_BITS),
        code.codeword_len(),
    );
    #[cfg(feature = "ligero")]
    tracing::info!(
        "ligero RS[n={}, k={}] δ={:.2}: {} of {} columns queried ({}-bit target, capacity conjecture 1−ρ)",
        LIGERO_RATE * cols,
        cols,
        1.0 - 1.0 / LIGERO_RATE as f64,
        num_queries,
        code.codeword_len(),
        TARGET_SECURITY_BITS,
    );
    let tensor = TensorPcs::new(code, mmcs, log_rows, num_queries);
    CodePcs::new(tensor)
}

/// A config for an AIR with a preprocessed trace: it carries two PCS instances,
/// one per stacked arity (main trace vs preprocessed trace).
pub struct EdDsaConfig {
    main_pcs: CodePcs,
    preprocessed_pcs: CodePcs,
}

impl MultiStarkConfig for EdDsaConfig {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = CodePcs;

    fn pcs(&self) -> &CodePcs {
        &self.main_pcs
    }

    fn preprocessed_pcs(&self) -> &CodePcs {
        &self.preprocessed_pcs
    }

    fn min_num_variables(&self) -> usize {
        BRAKEDOWN_FOLDING
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Witness<F> {
        L::new_witness(tables, BRAKEDOWN_FOLDING)
    }

    #[cfg(not(feature = "basefold"))]
    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_brakedown::BrakedownProverData<F, EF, MyMmcs, L>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }

    #[cfg(feature = "basefold")]
    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_basefold::BasefoldMlProverData<F, MyMmcs, L>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }
}

/// Build the EdDSA config: both traces are `2^log_height` rows tall; the main
/// trace is `main_width` columns and the preprocessed trace `preprocessed_width`.
pub fn config_for_eddsa(
    log_height: usize,
    main_width: usize,
    preprocessed_width: usize,
) -> EdDsaConfig {
    config_for_eddsa_with_security(log_height, main_width, preprocessed_width, TARGET_SECURITY_BITS)
}

/// As [`config_for_eddsa`], but with an explicit Brakedown soundness target. Production uses
/// [`TARGET_SECURITY_BITS`]; small traces (e.g. tests) may need a lower target, since secure
/// Brakedown at 100 bits requires enough columns to open `t` queries — a large `N`. Passing a
/// target the trace cannot support makes `build_pcs` panic rather than prove insecurely.
pub fn config_for_eddsa_with_security(
    log_height: usize,
    main_width: usize,
    preprocessed_width: usize,
    security_bits: usize,
) -> EdDsaConfig {
    EdDsaConfig {
        main_pcs: build_pcs(log_height, main_width, security_bits),
        preprocessed_pcs: build_pcs(log_height, preprocessed_width, security_bits),
    }
}
