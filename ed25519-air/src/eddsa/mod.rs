//! A single monolithic AIR proving ed25519 signature verification.
//!
//! One main trace and one preprocessed trace, no lookup or permutation arguments: every row type
//! (scalar-mul ladder step, subgroup check, final addition) shares the same columns, selected by
//! a preprocessed selector, and all cross-row wiring rides on carried columns.
//!
//! See `docs/EdDSA.md` for the row/column layout and the list of known gaps — notably,
//! the challenge scalar `k` enters as a public value and its relation to the message hash is
//! *not* proved here.

mod air;
mod columns;
mod generation;

pub use air::EdDsaAir;
pub use columns::{
    EDDSA_REGION_WIDTH, EdDsaLadderCols, EdDsaPreprocessedCols, EdDsaSharedCols,
    FIXED_LADDER_START, NUM_EDDSA_COLS, NUM_EDDSA_LADDER_COLS, NUM_EDDSA_PREPROCESSED_COLS,
    NUM_EDDSA_PUBLIC_VALUES, NUM_EDDSA_SHARED_COLS, SHARED_START, VAR_LADDER_START, add_cols,
    add_cols_mut, ladder_cols, ladder_cols_mut, sgc_cols, sgc_cols_mut, shared_cols,
    shared_cols_mut,
};
pub use generation::{
    add_row, bits_to_field, double_and_add, generate_preprocessed_trace,
    generate_preprocessed_trace_batch, generate_trace_rows,
    ladder_results, public_values, sgc_a_row, sgc_r_row, trace_height,
};
