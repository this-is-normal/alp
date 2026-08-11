mod air;
mod columns;
mod generation;

pub use air::EcAddAir;
pub use columns::{EcAddCols, NUM_EC_ADD_COLS};
pub use generation::{Point, add_points, generate_trace_rows};
