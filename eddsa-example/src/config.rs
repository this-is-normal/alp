//! Multi-stark backend for the EdDSA AIR, over the ed25519 base field and
//! committed with the Brakedown tensor PCS.
//!
//! This mirrors `ed25519-fib-example`'s config, but `EdDsaAir` declares a
//! **preprocessed trace** (the fixed-base doubling ladder + row selectors), so
//! the config must provide a *second* commitment scheme via
//! [`MultiStarkConfig::preprocessed_pcs`]. The preprocessed trace stacks a
//! different column count than the main trace (`8` vs `45`), giving it a
//! different stacked arity, so it needs its own Brakedown PCS instance sized for
//! that arity.

#[cfg(not(feature = "ligero"))]
use p3_brakedown::BrakedownCode as RealBrakedownCode;
#[cfg(feature = "ligero")]
use p3_brakedown::ReedSolomonCode;
use p3_brakedown::{BrakedownMlPcs, LinearCode, TensorPcs};
#[cfg(feature = "basefold")]
use p3_basefold::{BasefoldMlPcs, BasefoldPcs, RandomFoldableCode};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multi_stark::config::MultiStarkConfig;
use p3_sumcheck::layout::{Layout, PrefixProver, Table, Witness};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
#[cfg(feature = "count-ops")]
use p3_symmetric::{CryptographicPermutation, Permutation};
use p3_util::log2_ceil_usize;
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ---- field / hash / commitment stack (ed25519 + Poseidon2 + Merkle) ----

/// Base (trace) field.
pub type F = Ed25519;
/// Challenge field: the degree-1 self-extension of ed25519 (`|EF| ≈ 2^255`).
pub type EF = Ed25519;

/// Poseidon2 state width (ed25519 Poseidon2 supports width 3).
const WIDTH: usize = 3;
/// Sponge rate (challenger + hash); capacity = 1 (~255 bits).
const RATE: usize = 2;
/// Merkle digest / sponge output size in field elements.
const OUT: usize = 1;

#[cfg(not(feature = "count-ops"))]
type Perm = Poseidon2Ed25519<WIDTH>;
#[cfg(feature = "count-ops")]
type Perm = CountingPerm;

/// Counts every Poseidon2-over-Ed25519 permutation (a hash invocation). Locally the hash is over
/// Ed25519, so its muls are also counted in `MUL`; counting perms lets us subtract them to isolate
/// the PCS *arithmetic* muls (on vultr2/guest the hash is BabyBear, so this subtraction is moot).
#[cfg(feature = "count-ops")]
#[derive(Clone, Debug)]
pub struct CountingPerm(Poseidon2Ed25519<WIDTH>);
#[cfg(feature = "count-ops")]
impl Permutation<[Ed25519; WIDTH]> for CountingPerm {
    fn permute_mut(&self, input: &mut [Ed25519; WIDTH]) {
        p3_ed25519::counters::on_perm();
        self.0.permute_mut(input);
    }
}
#[cfg(feature = "count-ops")]
impl CryptographicPermutation<[Ed25519; WIDTH]> for CountingPerm {}
type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;

/// The general-`PrimeField` duplex challenger (ed25519 is not `PrimeField64`).
pub type MyChallenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;

type L = PrefixProver<F, EF>;

/// The tensor PCS's linear code. Default = the recursive **Brakedown** code
/// (linear-time). With `--features ligero` it becomes a **Reed–Solomon** code
/// (high distance → few queries, but superlinear Vandermonde encoding).
#[cfg(not(feature = "ligero"))]
pub type Code = RealBrakedownCode<F, EF>;
#[cfg(feature = "ligero")]
pub type Code = ReedSolomonCode<F>;

/// The multilinear PCS used by multi-stark. Default = the Brakedown tensor PCS; with
/// `--features basefold` = the ed25519-native Basefold PCS (polylog verifier).
#[cfg(not(feature = "basefold"))]
pub type CodePcs = BrakedownMlPcs<F, EF, Code, MyMmcs, L>;
#[cfg(feature = "basefold")]
pub type CodePcs = BasefoldMlPcs<F, MyMmcs, L>;

/// Fixed RNG seed for the (DEMO-grade) Poseidon2 round constants; seeding
/// deterministically makes the prover and verifier build the identical permutation.
const POSEIDON2_SEED: u64 = 0xED25519_C0FFEE;
/// Fixed seed for the Brakedown code's random expander matrices, so prover and
/// verifier build the identical linear code.
#[cfg(not(feature = "ligero"))]
const BRAKEDOWN_CODE_SEED: u64 = 0xB2A4E_D0;
/// Reed–Solomon inverse rate `n/k` (distance δ = 1 − 1/rate). Rate 2 → δ = 0.5.
#[cfg(feature = "ligero")]
const LIGERO_RATE: usize = 2;
/// Target security for the Ligero query count (bits).
#[cfg(feature = "ligero")]
const TARGET_SECURITY_BITS: usize = 100;
/// Column queries for the Brakedown path. DEMO-grade (functional, not secure):
/// at these sizes Brakedown's required query count (thousands) exceeds the
/// codeword, so no secure regime exists — see notes.
#[cfg(not(feature = "ligero"))]
const DEMO_QUERIES: usize = 32;

/// Number of column queries for a Reed–Solomon (Ligero) tensor test at
/// `security_bits` bits, under the **list-decoding-up-to-capacity** conjecture.
///
/// RS is MDS: coding rate `ρ = 1/rate`, relative distance `δ = 1 − ρ`. The
/// conjecture is that proximity can be tested up to the list-decoding *capacity*
/// radius `e = 1 − ρ = δ` (vs the proven unique-decoding `δ/2` or Johnson
/// `1 − √ρ`). A word `e`-far from the code then fails a random column check with
/// probability ≥ `e`, so `t` columns miss with probability `(1−e)^t = ρ^t`.
/// Setting `ρ^t ≤ 2^{−λ}` gives the clean
///
/// ```text
///   t = ⌈ λ / −log₂(1 − δ) ⌉ = ⌈ λ / log₂(rate) ⌉.
/// ```
///
/// This is the aggressive (conjectural) choice production STARKs use for concrete
/// parameters; it is NOT unconditionally proven.
#[cfg(feature = "ligero")]
fn ligero_num_queries(rate: usize, security_bits: usize) -> usize {
    let delta = 1.0 - 1.0 / rate as f64; // = 1 − ρ, list-decoding capacity radius
    let bits_per_query = -(1.0 - delta).log2(); // = log2(rate)
    (security_bits as f64 / bits_per_query).ceil() as usize
}
/// Padding floor. The tensor PCS commits at natural arity, so any positive arity is fine.
const BRAKEDOWN_FOLDING: usize = 1;

/// Build the fixed Poseidon2 permutation (ed25519 is not `PrimeField64`, so the
/// generic `new_from_rng` is used).
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

/// Calibration: ed25519 muls consumed by ONE Poseidon2 permutation (to subtract hashing from the
/// verify mul-count and isolate PCS arithmetic).
#[cfg(feature = "count-ops")]
pub fn muls_per_perm() -> u64 {
    use p3_field::PrimeCharacteristicRing;
    let p = perm();
    let mut st = [F::ZERO; WIDTH];
    p3_ed25519::counters::reset();
    p.permute_mut(&mut st);
    p3_ed25519::counters::snapshot().mul
}

/// Build one Brakedown PCS sized for a `2^log_height`-row, `width`-column trace.
///
/// The committed polynomial stacks every trace column, so the tensor scheme is
/// sized for the stacked variable count `log_height + ceil(log2(width))`, split
/// `log_rows = stacked / 2` row variables and the rest as columns.
#[cfg(feature = "basefold")]
fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    // Basefold commits the stacked polynomial (num_variables = log_height + ceil(log2 width)).
    let num_variables = log_height + log2_ceil_usize(width);
    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    // Rate 1/4 foldable code, ChaCha-seeded (cross-platform).
    let code = RandomFoldableCode::<F>::new(num_variables, 4, 0xBA5E_F01D);
    // 100-bit AGGRESSIVE security (rate 1/4, list-decoding radius δ≈0.215):
    // q = ⌈100 / −log₂(1−δ)⌉ = 288, clamped to the codeword length (= inv_rate·2^num_variables)
    // so tiny instances stay well-formed.
    let num_queries = core::cmp::min(288, 4usize << num_variables);
    BasefoldMlPcs::new(BasefoldPcs::new(code, mmcs, num_queries))
}

#[cfg(not(feature = "basefold"))]
fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    let stacked_num_variables = log_height + log2_ceil_usize(width);
    let log_rows = stacked_num_variables / 2;
    let cols = 1 << (stacked_num_variables - log_rows);

    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);

    #[cfg(not(feature = "ligero"))]
    let code = {
        // Recursive Brakedown code; seeded so prover and verifier build the identical code.
        let mut code_rng = SmallRng::seed_from_u64(BRAKEDOWN_CODE_SEED);
        RealBrakedownCode::from_message_len(cols, &mut code_rng)
    };
    #[cfg(feature = "ligero")]
    // Reed–Solomon `[n, k]` with k = cols (power of two) and n = rate·cols.
    let code = ReedSolomonCode::new(cols, LIGERO_RATE * cols);

    #[cfg(not(feature = "ligero"))]
    let num_queries = core::cmp::min(DEMO_QUERIES, code.codeword_len());
    #[cfg(feature = "ligero")]
    // Security-derived: t = ⌈λ / −log₂(1 − δ/2)⌉, clamped to the codeword length.
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

/// A Brakedown-PCS config for an AIR with a preprocessed trace: it carries two
/// PCS instances, one per stacked arity (main trace vs preprocessed trace).
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
    EdDsaConfig {
        main_pcs: build_pcs(log_height, main_width),
        preprocessed_pcs: build_pcs(log_height, preprocessed_width),
    }
}

/// A bare seeded challenger — identical on the prove and verify sides (Brakedown
/// injects no domain separator).
pub fn challenger() -> MyChallenger {
    MyChallenger::new(perm())
}
