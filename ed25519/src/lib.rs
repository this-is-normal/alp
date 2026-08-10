//! The ed25519 base field, defined as `F_p` where `p = 2^255 - 19`.
//!
//! This is a 255-bit prime field implemented in the same style as [`p3_bn254`], using a
//! 4-limb little-endian Montgomery representation with respect to `R = 2^256`.
//!
//! Unlike bn254, this field is **not** two-adic-friendly (`p - 1 = 2^2 * (2^253 - 5)` has
//! 2-adicity only 2), so it does **not** implement `TwoAdicField`. It is also not
//! `PrimeField64`/`PrimeField32`.
#![no_std]

extern crate alloc;

mod ed25519;
mod helpers;
mod poseidon2;

pub use ed25519::*;
pub use poseidon2::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, ED25519_S_BOX_DEGREE,
    Poseidon2Ed25519, Poseidon2ExternalLayerEd25519, Poseidon2InternalLayerEd25519,
};
