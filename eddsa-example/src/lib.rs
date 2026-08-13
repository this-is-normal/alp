//! End-to-end ed25519 signature-verification STARK: the monolithic
//! [`p3_ed25519_air::eddsa::EdDsaAir`] proved over the **native ed25519 base
//! field** with the field-agnostic Brakedown tensor PCS on the multilinear
//! `multi-stark` path.
//!
//! * [`config`] — the ed25519 + Poseidon2 + Merkle + Brakedown backend, with the
//!   second (preprocessed) commitment scheme `EdDsaAir` requires.
//! * [`curve`] — real ed25519 constants and a genuine verifying instance.

pub mod config;
pub mod curve;
