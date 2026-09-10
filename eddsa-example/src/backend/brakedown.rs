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
/// Target soundness for the column-query count.
const SECURITY_BITS: usize = 100;
/// `false` = provable δ/3 radius (Ligero/Roth–Zémor); `true` = conjectured δ/2
/// (Diamond–Posen Conjecture-1).
const CONJECTURED: bool = false;

/// Build one tensor PCS sized for a `2^log_height`-row, `width`-column trace.
///
/// The committed polynomial stacks every trace column, so the scheme is sized
/// for the stacked variable count `log_height + ceil(log2(width))`, split by the
/// paper's optimal row/column rule.
pub fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    let stacked = log_height + log2_ceil_usize(width);

    // Small-proof preset (Haböck SdigCode6, δ=0.070, weights 9/12) with a
    // distance-derived query count.
    let preset = p3_brakedown::params::small_proof();
    let t = preset.num_queries(SECURITY_BITS, CONJECTURED);
    let log_rows = p3_brakedown::params::optimal_log_rows(stacked, t, &preset.params);
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
