mod air;
mod columns;
mod generation;

pub use air::EcScalarMulAir;
pub use columns::{EcScalarMulCols, NUM_EC_SCALARMUL_COLS};
pub use generation::generate_trace_rows;
