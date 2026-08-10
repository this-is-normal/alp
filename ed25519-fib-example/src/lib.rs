//! End-to-end **native** multilinear STARK over the ed25519 base field
//! (`p = 2^255 - 19`), committed with the field-agnostic Brakedown tensor PCS.
//!
//! This mirrors `multi-fib-example`'s Brakedown path, but instantiated entirely
//! over the ed25519 field, which is **non-two-adic** and neither `PrimeField64`
//! nor `PrimeField32`. It demonstrates that the multi-stark + Brakedown stack is
//! genuinely field-agnostic:
//!
//! * `F = Ed25519` is both the trace field **and** the challenge field
//!   (`EF = Ed25519`, a degree-1 self-extension — soundness is fine since
//!   `|EF| ≈ 2^255`).
//! * The Fiat–Shamir transcript uses [`DuplexChallengerPF`], the general
//!   `PrimeField` challenger wrapper (the bare `DuplexChallenger`'s bit/grinding
//!   impls require `PrimeField64`).
//! * Merkle commitments use a width-3 Poseidon2 over ed25519.

pub mod config;
pub mod fib_air;
