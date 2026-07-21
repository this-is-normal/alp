# Project Structure

A map of the crates and folders in this repository (a fork of Plonky3 tailored
to EdDSA verification). Descriptions are drawn from each crate's `Cargo.toml`
and `README.md`.

The workspace is a toolkit of primitives for building polynomial IOPs and
STARK-based proof systems. Components are grouped below by role.

## Finite fields

- **`field/`** — Modular framework for finite fields: generic binomial
  extensions and SIMD-packed field arithmetic. The core traits everything else
  builds on.
- **`field-testing/`** — Shared testing utilities for finite-field arithmetic
  and field-based operations.
- **`monty-31/`** — Generic prime field `F_p` for `2^30 < p < 2^31` using
  Montgomery arithmetic; the shared backbone for BabyBear and KoalaBear.
- **`baby-bear/`** — The BabyBear prime field `F_p`, `p = 2^31 - 2^27 + 1`.
- **`koala-bear/`** — The KoalaBear prime field `F_p`, `p = 2^31 - 2^24 + 1`.
- **`mersenne-31/`** — The Mersenne31 prime field `F_p`, `p = 2^31 - 1`.
- **`goldilocks/`** — The Goldilocks prime field `F_p`, `p = 2^64 - 2^32 + 1`.
- **`bn254/`** — Wrapper around `ark_bn254::Fr` so it can be used as a Plonky3
  field.

## Linear algebra and transforms

- **`matrix/`** — Matrix library with common operations and traits over
  finite-field elements.
- **`dft/`** — Discrete Fourier transform (DFT) implementations for finite
  fields.
- **`mds/`** — Fast multiplication by MDS matrices (used inside permutations).
- **`interpolation/`** — Polynomial interpolation utilities (currently a
  placeholder crate).
- **`multilinear-util/`** — Utilities for multilinear polynomial operations.

## Symmetric primitives and permutations

- **`symmetric/`** — Modular framework for symmetric cryptographic primitives
  (sponges, compression functions, hashers).
- **`poseidon1/`**, **`poseidon2/`** — The Poseidon1 and Poseidon2
  permutations for hashing in ZK applications.
- **`monolith/`** — The Monolith permutation family over Mersenne31 and
  Goldilocks.
- **`rescue/`** — The Rescue hash function and its variants.
- **`keccak/`** — Optimized Keccak-f permutation and Keccak256 hash.
- **`blake3/`** — Wrapper around the standard BLAKE3 hash for use in Plonky3.
- **`sha256/`** — Plonky3 hash-trait implementations for SHA2-256.

## AIRs (Algebraic Intermediate Representations)

- **`air/`** — The AIR interface and supporting utilities that STARK provers
  target.
- **`poseidon1-air/`**, **`poseidon2-air/`** — AIRs for the Poseidon1 /
  Poseidon2 permutations.
- **`monolith-air/`** — AIR for the Monolith permutation family.
- **`keccak-air/`** — AIR for the Keccak-f permutation (field size `2^16`–`2^32`).
- **`blake3-air/`** — AIR for the Blake3 permutation (field size `2^20`–`2^32`).

## Commitment schemes and low-degree testing

- **`commit/`** — Framework for cryptographic commitment schemes, including
  non-hiding variants.
- **`merkle-tree/`** — Merkle-tree commitment scheme.
- **`fri/`** — The FRI low-degree test (LDT) protocol.
- **`whir/`** — The WHIR polynomial commitment scheme.
- **`zk-codes/`** — Zero-knowledge codes, including Reed–Solomon ZK encoding.

## Proof systems and protocol engines

- **`uni-stark/`** — Minimal univariate STARK framework: core proof generation
  and verification.
- **`batch-stark/`** — Batched STARK wrapper over `uni-stark` that reuses FRI
  openings across instances.
- **`multi-stark/`** — Multilinear (SuperSpartan-style) STARK prover for AIRs.
- **`circle/`** — Circle STARK proof system built around the unit circle of a
  finite field.
- **`sumcheck/`** — Sumcheck protocol engine (product polynomials, SVO, HVZK
  masks).
- **`lookup/`** — LogUp: logarithmic lookup arguments and multiset-equality
  checks with multiplicities.
- **`challenger/`** — Fiat–Shamir transcript / challenger framework for
  deriving random challenges.

## Support, tooling, and analysis

- **`util/`** — Low-level utilities (bit manipulation, array transformations).
- **`maybe-rayon/`** — Feature-gated wrapper around `rayon` for optional
  parallelism.
- **`examples/`** — Example implementations demonstrating how the components
  fit together.
- **`security/`** — Soundness analysis for STARK protocols: per-error-source
  modules (FRI, DEEP, grinding, proximity, …) composed at the protocol call
  site.
- **`audits/`** — External audit reports (e.g. the Least Authority / Polygon
  Plonky3 report).
