mod air;
mod columns;
mod generation;

pub use air::{EcScalarMulAir, NUM_PUBLIC_VALUES};
pub use columns::{EcScalarMulCols, NUM_EC_SCALARMUL_COLS};
pub use generation::{ScalarMulTrace, generate_trace_rows};
