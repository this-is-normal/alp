//! **Ligero**: the same tensor PCS over a Reed–Solomon code. High distance means
//! few queries and small proofs, at the cost of superlinear Vandermonde
//! encoding.

use p3_brakedown::{BrakedownMlPcs, LinearCode, ReedSolomonCode, TensorPcs};
use p3_util::log2_ceil_usize;

use crate::hash::{EF, F, L, MyMmcs, build_mmcs};

pub type Code = ReedSolomonCode<F>;
/// Human-readable backend name, for measurement output.
pub const NAME: &str = "Ligero";

pub type CodePcs = BrakedownMlPcs<F, EF, Code, MyMmcs, L>;

/// Inverse rate `n/k`; distance `δ = 1 − 1/rate`, so rate 2 gives `δ = 0.5`.
const RATE: usize = 2;
/// Target security for the query count (bits).
const SECURITY_BITS: usize = 100;

/// Query count under the **list-decoding-up-to-capacity** conjecture:
/// `t = ⌈λ / log₂(rate)⌉`. This is the aggressive choice production STARKs use
/// for concrete parameters; it is NOT unconditionally proven.
fn num_queries(rate: usize, security_bits: usize) -> usize {
    let delta = 1.0 - 1.0 / rate as f64; // = 1 − ρ, the capacity radius
    (security_bits as f64 / -(1.0 - delta).log2()).ceil() as usize
}

/// Build one RS tensor PCS sized for a `2^log_height`-row, `width`-column trace,
/// with a square row/column split.
pub fn build_pcs(log_height: usize, width: usize) -> CodePcs {
    let stacked = log_height + log2_ceil_usize(width);
    let log_rows = stacked / 2;
    let cols = 1 << (stacked - log_rows);

    let code = ReedSolomonCode::new(cols, RATE * cols);
    let t = core::cmp::min(num_queries(RATE, SECURITY_BITS), code.codeword_len());
    tracing::info!(
        "ligero RS[n={}, k={}] δ={:.2}: {} of {} columns queried ({}-bit target, capacity conjecture 1−ρ)",
        RATE * cols,
        cols,
        1.0 - 1.0 / RATE as f64,
        t,
        code.codeword_len(),
        SECURITY_BITS,
    );
    CodePcs::new(TensorPcs::new(code, build_mmcs(), log_rows, t))
}
