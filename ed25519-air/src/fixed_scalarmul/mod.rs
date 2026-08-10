mod air;
mod columns;
mod generation;

pub use air::EcFixedScalarMulAir;
pub use columns::{
    EcFixedScalarMulCols, EcFixedScalarMulPreprocessedCols, NUM_EC_FIXED_SCALARMUL_COLS,
    NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS,
};
pub use generation::{generate_preprocessed_trace, generate_trace_rows};
