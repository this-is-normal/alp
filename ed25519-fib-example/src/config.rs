//! The multi-stark backend wiring over the ed25519 base field, committed with
//! the Brakedown tensor PCS.
//!
//! This mirrors `multi-fib-example`'s `MultiFibConfigBrakedown`, but with:
//! * `F = EF = Ed25519` (a non-two-adic, non-`PrimeField64` field used as both
//!   the trace field and the degree-1 self-extension challenge field),
//! * a width-3 Poseidon2 over ed25519 for the Merkle hash/compression,
//! * the general-`PrimeField` [`DuplexChallengerPF`] transcript wrapper.

use p3_brakedown::{BrakedownCode as RealBrakedownCode, BrakedownMlPcs, LinearCode, TensorPcs};
use p3_challenger::DuplexChallengerPF;
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Ed25519, Poseidon2Ed25519,
};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_multi_stark::config::MultiStarkConfig;
use p3_sumcheck::layout::{Layout, PrefixProver, Table, Witness};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_util::log2_ceil_usize;
use rand::SeedableRng;
use rand::rngs::SmallRng;

// ---- field / hash / commitment stack (ed25519 + Poseidon2 + Merkle) ----

/// Base (trace) field.
pub type F = Ed25519;
/// Challenge field: the degree-1 self-extension of ed25519. `|EF| ≈ 2^255` gives
/// ample soundness without needing a genuine field extension.
pub type EF = Ed25519;

/// Poseidon2 state width. bn254/ed25519 Poseidon2 only supports width 3.
const WIDTH: usize = 3;
/// Sponge rate (challenger + hash). `RATE < WIDTH`, capacity = 1 (~255 bits).
const RATE: usize = 2;
/// Merkle digest / sponge output size in field elements. `2*OUT <= WIDTH` (binary
/// tree via `TruncatedPermutation`), and `OUT <= RATE`. One 255-bit element gives
/// ~127-bit collision resistance.
const OUT: usize = 1;

type Perm = Poseidon2Ed25519<WIDTH>;

type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;

/// The general-`PrimeField` duplex challenger (ed25519 is not `PrimeField64`).
pub type MyChallenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;

type L = PrefixProver<F, EF>;

/// The real recursive Brakedown linear-time code the tensor PCS commits through.
pub type BrakedownCode = RealBrakedownCode<F, EF>;
/// The standalone tensor commitment wrapped by the multilinear adapter.
pub type BrakedownTensor = TensorPcs<F, EF, BrakedownCode, MyMmcs>;
/// The Brakedown-backed multilinear PCS used by multi-stark.
pub type BrakedownPcs = BrakedownMlPcs<F, EF, BrakedownCode, MyMmcs, L>;

/// Fixed RNG seed for the Poseidon2 round constants (DEMO-grade; see
/// `p3_ed25519::poseidon2`). Seeding deterministically makes the prover and
/// verifier build the identical permutation.
const POSEIDON2_SEED: u64 = 0xED25519_C0FFEE;

/// Fixed seed for the Brakedown code's random `A`/`B` expander matrices, so the
/// prover and verifier configs build the **identical** linear code.
const BRAKEDOWN_CODE_SEED: u64 = 0xB2A4E_D0;

/// Preprocessing/padding floor. The tensor PCS commits the stacked polynomial at
/// its natural arity, so any positive arity is fine.
const BRAKEDOWN_FOLDING: usize = 1;

/// Build the fixed Poseidon2 permutation. Uses the generic `new_from_rng`
/// (ed25519 is not `PrimeField64`, so `new_from_rng_128` is unavailable).
fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(POSEIDON2_SEED);
    Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    )
}

/// A Brakedown-PCS-backed multilinear AIR configuration over ed25519.
pub struct MultiFibConfigBrakedown {
    pcs: BrakedownPcs,
}

impl MultiStarkConfig for MultiFibConfigBrakedown {
    type Val = F;
    type Challenge = EF;
    type Challenger = MyChallenger;
    type Pcs = BrakedownPcs;

    fn pcs(&self) -> &BrakedownPcs {
        &self.pcs
    }

    fn min_num_variables(&self) -> usize {
        BRAKEDOWN_FOLDING
    }

    fn build_witness(&self, tables: Vec<Table<F>>) -> Witness<F> {
        L::new_witness(tables, BRAKEDOWN_FOLDING)
    }

    fn committed_table<'a>(
        &self,
        prover_data: &'a p3_brakedown::BrakedownProverData<F, EF, MyMmcs, L>,
        table_index: usize,
    ) -> &'a Table<F> {
        prover_data.table(table_index)
    }
}

/// Build a Brakedown-PCS configuration sized for a single trace of
/// `2^log_height` rows and `width` columns.
///
/// The committed polynomial stacks every trace column, so the tensor scheme is
/// sized for the stacked variable count `log_height + ceil(log2(width))`, split
/// into `log_rows = stacked / 2` row variables and the rest as columns.
pub fn config_for_brakedown(log_height: usize, width: usize) -> MultiFibConfigBrakedown {
    let stacked_num_variables = log_height + log2_ceil_usize(width);
    let log_rows = stacked_num_variables / 2;
    let cols = 1 << (stacked_num_variables - log_rows);

    let mmcs = MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), 0);
    // Real recursive Brakedown code over the (power-of-two) message length `cols`.
    let mut code_rng = SmallRng::seed_from_u64(BRAKEDOWN_CODE_SEED);
    let code = BrakedownCode::from_message_len(cols, &mut code_rng);
    // Brakedown's codeword length is not a power of two; clamp the query count to it.
    let num_queries = core::cmp::min(32, code.codeword_len());
    let tensor = TensorPcs::new(code, mmcs, log_rows, num_queries);
    MultiFibConfigBrakedown {
        pcs: BrakedownPcs::new(tensor),
    }
}

/// The Brakedown PCS injects no domain separator, so a bare seeded challenger
/// suffices — identical on the prove and verify sides.
pub fn challenger_brakedown() -> MyChallenger {
    MyChallenger::new(perm())
}
