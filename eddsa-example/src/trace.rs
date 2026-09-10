//! Building a batched EdDSA trace: `2^log_instances` signatures stacked
//! vertically into one table.
//!
//! `EdDsaAir` supports this natively — its preprocessed trace restarts the
//! selectors (`is_instance_first` / `is_instance_last`) at every instance
//! boundary, and the per-instance constraints are gated on those rather than on
//! the global first row — but the AIR crate ships only the single-signature
//! main-trace generator, so the stacking lives here.
//!
//! # Scope
//!
//! Every stacked instance is *fully constrained* as a signature verification.
//! Only instance 0 is bound to the proof's public values (the AIR reads them
//! `when_first_row`), so a batched proof says "I know 2^k valid signatures, the
//! first of which is this one" rather than exposing all 2^k public keys. That is
//! the right shape for a cost comparison, and it is the same boundary the
//! single-signature path already has (`k = H(...)` is not proved either).

use p3_ed25519_air::eddsa::{NUM_EDDSA_COLS, generate_trace_rows};
use p3_matrix::dense::RowMajorMatrix;

use crate::curve::{Instance, base_b, d, subgroup_order};
use crate::hash::F;

/// Stack one main trace per instance. The result is
/// `instances.len() * trace_height(num_bits)` rows tall.
///
/// # Panics
///
/// Panics if `instances` is empty, or if the per-instance traces disagree in
/// height (they cannot: all share `num_bits`).
#[must_use]
pub fn generate_batch_trace(instances: &[Instance]) -> RowMajorMatrix<F> {
    assert!(!instances.is_empty(), "need at least one signature");
    let (d, base_b, l) = (d(), base_b(), subgroup_order());

    let mut values = Vec::new();
    for instance in instances {
        let trace = generate_trace_rows(
            d,
            base_b,
            &l,
            instance.a,
            instance.r,
            &instance.s_bits,
            &instance.k_bits,
        );
        values.extend_from_slice(&trace.values);
    }
    RowMajorMatrix::new(values, NUM_EDDSA_COLS)
}
