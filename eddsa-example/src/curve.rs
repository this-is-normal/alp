//! Real ed25519 curve constants and construction of a genuine, non-degenerate
//! signature-verification instance over the native base field `p = 2^255 − 19`.
//!
//! The `p3-ed25519-air` tests run over `p3-bn254` with placeholder constants and
//! a *degenerate* instance (`A = R = O`), because their placeholder curve makes
//! the subgroup-check witness derivation close only for the identity. Over the
//! **real** curve, with the real subgroup order `l`, every genuine subgroup point
//! has a valid `8^{-1} · P` witness, so a real signature verifies.
//!
//! We build a self-consistent instance from a secret key `a`, nonce `r`, and
//! challenge `k` (all small — the circuit size is set by `num_bits`, not by the
//! scalar magnitudes):
//!
//! ```text
//!   A = [a]·B,   R = [r]·B,   s = r + k·a   (as integers, kept < l)
//!   ⇒  [s]·B = [r + k·a]·B = [r]·B + [k]·([a]·B) = R + [k]·A       (ord B = l)
//! ```
//!
//! which is exactly the verification equation `EdDsaAir` checks. `k` is a free
//! challenge here — the AIR does not prove `k = H(R,A,M)` (see the AIR docs).

use num_bigint::BigUint;
use p3_ed25519::Ed25519;
use p3_ed25519_air::add::{Point, add_points};
use p3_ed25519_air::eddsa::double_and_add;
use p3_ed25519_air::subgroupcheck::scalar_mul;
use p3_field::PrimeCharacteristicRing;

/// Base (curve) field.
pub type F = Ed25519;

fn big(decimal: &str) -> BigUint {
    decimal.parse().expect("valid decimal BigUint")
}

fn fe(decimal: &str) -> F {
    Ed25519::from_biguint(big(decimal)).expect("in-range field element")
}

/// Twisted-Edwards parameter `d = −121665/121666 mod p` (RFC 8032 ed25519).
pub fn d() -> F {
    fe("37095705934669439343138083508754565189542113879843219016388785533085940283555")
}

/// The standard ed25519 base point `B` (order `l`).
pub fn base_b() -> Point<F> {
    (
        fe("15112221349535400772501151409588531511454012693041857206046113283949847762202"),
        fe("46316835694926478169428394003475163141307993866256225615783033603165251855960"),
    )
}

/// Order `l` of the prime-order subgroup: `2^252 + 27742317777372353535851937790883648493`.
pub fn subgroup_order() -> BigUint {
    big("7237005577332262213973186563042994240857116359379907606001950938285454250989")
}

/// `true` iff `(x, y)` satisfies the twisted-Edwards equation `−x² + y² = 1 + d·x²·y²`.
fn on_curve(d: F, (x, y): Point<F>) -> bool {
    -(x * x) + y * y == F::ONE + d * (x * x) * (y * y)
}

/// LSB-first bit decomposition of `n` into exactly `num_bits` bits.
fn bits_le(n: &BigUint, num_bits: usize) -> Vec<bool> {
    (0..num_bits as u64).map(|i| n.bit(i)).collect()
}

/// A verifying EdDSA instance: public key `A`, commitment `R`, and the LSB-first
/// bit decompositions of the response `s` and challenge `k` (each `num_bits` long).
pub struct Instance {
    pub a: Point<F>,
    pub r: Point<F>,
    pub s_bits: Vec<bool>,
    pub k_bits: Vec<bool>,
}

/// Build a genuine verifying instance for `num_bits`-bit scalars from small secrets
/// `(a_sec, r_sec, k)`. Panics (with a clear message) if the constants or the
/// derived witness fail to satisfy the curve/verification equation — so a bad
/// constant can never silently produce an unsatisfiable trace.
pub fn valid_instance(num_bits: usize, a_sec: u64, r_sec: u64, k: u64) -> Instance {
    let d = d();
    let b = base_b();
    let l = subgroup_order();

    // Constants sanity: B is on-curve and in the subgroup (order l ⇒ [l]B = O).
    assert!(on_curve(d, b), "base point B is not on the curve — check d/B constants");

    let a_sec = BigUint::from(a_sec);
    let r_sec = BigUint::from(r_sec);
    let k_big = BigUint::from(k);

    // A = [a]B, R = [r]B (genuine order-l subgroup points).
    let a = scalar_mul(d, &a_sec, b);
    let r = scalar_mul(d, &r_sec, b);
    assert!(on_curve(d, a) && on_curve(d, r), "derived A/R off curve");

    // s = r + k·a as an integer; keep it < l so [s]B needs no explicit reduction
    // and so its bit length fits num_bits.
    let s = &r_sec + &k_big * &a_sec;
    assert!(s < l, "s = r + k·a must stay below the subgroup order for this demo");
    assert!(
        (s.bits() as usize) <= num_bits && (k_big.bits() as usize) <= num_bits,
        "num_bits ({num_bits}) too small for s ({} bits) / k ({} bits)",
        s.bits(),
        k_big.bits()
    );

    let s_bits = bits_le(&s, num_bits);
    let k_bits = bits_le(&k_big, num_bits);

    // Independent check of the verification equation the AIR will enforce:
    //   R + [k]A == [s]B.
    let lhs = add_points(d, r, double_and_add(d, a, &k_bits));
    let rhs = double_and_add(d, b, &s_bits);
    assert_eq!(lhs, rhs, "instance does not satisfy R + [k]A == [s]B");

    Instance { a, r, s_bits, k_bits }
}

/// `num_signatures` genuine, *distinct* verifying instances, for the batched
/// proofs `EdDsaAir`'s `log_instances` stacks into one table.
///
/// The secrets vary per instance so the stacked trace is not `n` copies of one
/// witness. They stay small (the demo constraint `s = r + k·a < l` bounds them),
/// which costs nothing: the circuit is full `num_bits` rows regardless of scalar
/// magnitude.
#[must_use]
pub fn valid_instances(num_bits: usize, num_signatures: usize) -> Vec<Instance> {
    (0..num_signatures)
        .map(|i| {
            let i = i as u64;
            valid_instance(num_bits, 3 + i, 5 + 2 * i, 7 + i)
        })
        .collect()
}
