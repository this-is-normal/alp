//! End-to-end ed25519 signature-verification STARK: the monolithic
//! [`p3_ed25519_air::eddsa::EdDsaAir`] proved over the **native ed25519 base
//! field** on the multilinear `multi-stark` path, with a field-agnostic
//! commitment scheme chosen by cargo feature.
//!
//! * [`config`] — the ed25519 + Poseidon2 + Merkle + Brakedown backend, with the
//!   second (preprocessed) commitment scheme `EdDsaAir` requires.
//! * [`hash`] — the field / Poseidon2 / Merkle / transcript stack.
//! * [`backend`] — one commitment scheme per cargo feature (Brakedown, Ligero,
//!   BaseFold, Quasar).
//! * [`curve`] — real ed25519 constants and a genuine verifying instance.

pub mod backend;
pub mod config;
pub mod curve;
pub mod trace;
pub mod hash;
