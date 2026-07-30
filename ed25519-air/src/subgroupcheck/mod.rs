mod air;
mod columns;
mod generation;

pub use air::EcSubgroupCheckAir;
pub use columns::{EcSubgroupCheckCols, NUM_EC_SUBGROUPCHECK_COLS};
pub use generation::{
    generate_trace_rows, generate_trace_rows_for_targets, scalar_mul, witness_for_target,
};
