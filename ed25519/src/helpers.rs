use alloc::vec::Vec;

use num_bigint::BigUint;

use crate::{ED25519_MONTY_MU_64, ED25519_PRIME};

/// Const lexicographic comparison: returns true if a < b (little-endian limbs).
#[inline]
pub(crate) const fn const_lt(a: [u64; 4], b: [u64; 4]) -> bool {
    // Compare from most significant limb to least significant
    if a[3] != b[3] {
        return a[3] < b[3];
    }
    if a[2] != b[2] {
        return a[2] < b[2];
    }
    if a[1] != b[1] {
        return a[1] < b[1];
    }
    a[0] < b[0]
}

/// Convert a fixed-size array of u64s (little-endian) to a BigUint.
#[inline]
pub(crate) fn to_biguint<const N: usize>(value: [u64; N]) -> BigUint {
    let bytes: Vec<u8> = value.iter().flat_map(|x| x.to_le_bytes()).collect();
    BigUint::from_bytes_le(&bytes)
}

/// Basically copied the implementation here: <https://doc.rust-lang.org/std/primitive.u32.html#method.carrying_add>
///
/// Once this moves to standard rust (currently nightly) we can use that directly.
/// Tracking Issue is here: <https://github.com/rust-lang/rust/issues/85532>
#[inline]
const fn carrying_add(lhs: u64, rhs: u64, carry: bool) -> (u64, bool) {
    let (a, c1) = lhs.overflowing_add(rhs);
    let (b, c2) = a.overflowing_add(carry as u64);
    (b, c1 | c2)
}

/// Compute `lhs + rhs`, returning a bool if overflow occurred.
#[inline]
pub(crate) const fn wrapping_add(lhs: [u64; 4], rhs: [u64; 4]) -> ([u64; 4], bool) {
    let mut carry = false;
    let mut output = [0; 4];

    let mut i = 0;
    while i < 4 {
        (output[i], carry) = carrying_add(lhs[i], rhs[i], carry);
        i += 1;
    }

    (output, carry)
}

/// Basically copied the implementation here: <https://doc.rust-lang.org/std/primitive.u32.html#method.borrowing_sub>
///
/// Once this moves to standard rust (currently nightly) we can use that directly.
/// Tracking Issue is here: <https://github.com/rust-lang/rust/issues/85532>.
#[inline]
const fn borrowing_sub(lhs: u64, rhs: u64, borrow: bool) -> (u64, bool) {
    let (a, c1) = lhs.overflowing_sub(rhs);
    let (b, c2) = a.overflowing_sub(borrow as u64);
    (b, c1 | c2)
}

/// Compute `lhs - rhs`, returning a bool if underflow occurred.
#[inline]
pub(crate) const fn wrapping_sub(lhs: [u64; 4], rhs: [u64; 4]) -> ([u64; 4], bool) {
    let mut borrow = false;
    let mut output = [0; 4];

    // use manual `while` loop to enable `const`
    let mut i = 0;
    while i < 4 {
        (output[i], borrow) = borrowing_sub(lhs[i], rhs[i], borrow);
        i += 1;
    }

    (output, borrow)
}

/// Compute a * b with a in the range 0..2^256 and b in the range 0..2^64.
///
/// Returns the lowest output limb and the remaining limbs in a 4-limb array.
#[inline]
pub(crate) const fn mul_small(lhs: [u64; 4], rhs: u64) -> (u64, [u64; 4]) {
    let mut output = [0u64; 4];
    let mut acc;

    // Process the first limb separately to get the lowest output limb.
    acc = (lhs[0] as u128) * (rhs as u128);
    let out_0 = acc as u64;

    // acc < 2^64
    acc >>= 64;

    // Process the remaining limbs.
    // use manual `while` loop to enable `const`
    let mut i = 1;
    while i < 4 {
        // Product of u64's < 2^128 - 2^64 so this addition will not overflow.
        acc += (lhs[i] as u128) * (rhs as u128);
        output[i - 1] = acc as u64;

        // acc < 2^64
        acc >>= 64;

        i += 1;
    }
    output[3] = acc as u64;

    (out_0, output)
}

/// Compute a * b + c with a, c in the range 0..2^256 and b in the range 0..2^64.
///
/// Returns the lowest output limb and the remaining limbs in a 4-limb array.
#[inline]
pub(crate) const fn mul_small_and_acc(lhs: [u64; 4], rhs: u64, add: [u64; 4]) -> (u64, [u64; 4]) {
    let mut output = [0u64; 4];
    let mut acc;

    // Process the first limb separately to get the lowest output limb.
    acc = (lhs[0] as u128) * (rhs as u128) + (add[0] as u128);
    let out_0 = acc as u64;

    // acc < 2^64
    acc >>= 64;

    // Process the remaining limbs.
    // use manual `while` loop to enable `const`
    let mut i = 1;
    while i < 4 {
        // Product of u64's < 2^128 - 2^64 so this addition will not overflow.
        acc += (lhs[i] as u128) * (rhs as u128) + (add[i] as u128);
        output[i - 1] = acc as u64;

        // acc < 2^64
        acc >>= 64;

        i += 1;
    }
    output[3] = acc as u64;

    (out_0, output)
}

// Interleaved Montgomery multiplication:
//
// When working with Big-Nums where the base multiplication is expensive, we
// use a variant of montgomery multiplication which is more efficient. The idea
// is to interleave the multiplication and reduction steps which lets us
// avoid the need for Big-Num x Big-Num multiplications.
//
// Let P be our prime and `mu = P^{-1} mod 2^64`.
// The Interleaved Montgomery reduction (IMR) algorithm takes as input a 320-bit number `x`
// and returns a 256-bit number equal to `2^{-64}x mod P`.
//
// 1. Define `t = x * mu mod 2^64`.
// 2. Define `u = t * P`.
// 3. Define `sub = (x - u) / 2^{64}`.
// 4. If `sub < 0`, return `sub + P` else return `sub`.
//
// The division in step 3 is exact as `u mod 2^64 = t * P mod 2^64 = x * mu * P mod 2^64 = x mod 2^64`.
// Additionally note that the output is `< P` if `sub < 0` and otherwise is `< x/2^{64}`. Hence if we
// assume that `x < 2^{64}P`, then the output is always `< P`.
//
// This works for any odd prime `P < 2^255` (in particular `p = 2^255 - 19`), following the exact
// same induction as bn254: assuming `lhs < P` the running result stays `< P` after each of the four
// reduction steps.

/// The incoming number is split into the low 64-bit limb `acc0` (treated separately) and the
/// remaining four 64-bit limbs `acc`.
#[inline]
const fn interleaved_monty_reduction(acc0: u64, acc: [u64; 4]) -> [u64; 4] {
    let t = acc0.wrapping_mul(ED25519_MONTY_MU_64);
    let (_, u) = mul_small(ED25519_PRIME, t);

    let (sub, under) = wrapping_sub(acc, u);
    if under {
        let (sub_corr, _) = wrapping_add(sub, ED25519_PRIME);
        sub_corr
    } else {
        sub
    }
}

/// Montgomery multiplication and reduction algorithm for the ed25519 base field.
///
/// The algorithm assumes that `lhs < P` but puts no constraints on `rhs`.
///
/// The output is a 4-limb array representing the result of `lhs * rhs * 2^{-256} mod P`
/// guaranteed to be in the range `[0, P)`.
#[inline]
pub(crate) const fn monty_mul(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
    // We need to ensure that `lhs < P` otherwise it's possible for the
    // algorithm to fail and produce a value which is too large.
    debug_assert!(const_lt(lhs, ED25519_PRIME));

    // Our accumulator starts at 0 so we start with mul_small
    let (acc0, acc) = mul_small(lhs, rhs[0]);
    let res0 = interleaved_monty_reduction(acc0, acc);

    // Then we repeat the above process for the remaining rhs limbs
    // including the previous result as an accumulator.
    let (acc0, acc) = mul_small_and_acc(lhs, rhs[1], res0);
    let res1 = interleaved_monty_reduction(acc0, acc);
    let (acc0, acc) = mul_small_and_acc(lhs, rhs[2], res1);
    let res2 = interleaved_monty_reduction(acc0, acc);
    let (acc0, acc) = mul_small_and_acc(lhs, rhs[3], res2);
    interleaved_monty_reduction(acc0, acc)
}

/// The ed25519 prime `p = 2^255 - 19` represented as a little-endian array of 2-u128s.
const ED25519_PRIME_U128: [u128; 2] = [
    0xffffffffffffffffffffffffffffffed,
    0x7fffffffffffffffffffffffffffffff,
];

/// Efficiently halve an ed25519 field element.
#[inline]
pub(crate) const fn halve_ed25519(mut input: [u64; 4]) -> [u64; 4] {
    // Seems to be a little faster to convert into u128s.
    // It's essentially identical under the hood so this is
    // likely just helping the compiler generate simpler assembly somehow.
    let mut input0_u128 = (input[1] as u128) << 64 | (input[0] as u128);
    let mut input1_u128 = (input[3] as u128) << 64 | (input[2] as u128);
    let carry;

    if input0_u128 & 1 == 1 {
        // If the element is odd, we add P (which is odd), making it even.
        (input0_u128, carry) = input0_u128.overflowing_add(ED25519_PRIME_U128[0]);
        // Can ignore overflow here as the sum is < 2^256.
        input1_u128 = input1_u128.wrapping_add(ED25519_PRIME_U128[1]);
        input1_u128 = input1_u128.wrapping_add(carry as u128);
    }

    // At this point the element is guaranteed to be even so we just
    // need to shift down by 1.
    let carry_bit = (input1_u128 << 63) as u64;

    input[0] = (input0_u128 >> 1) as u64;
    input[1] = (input0_u128 >> 65) as u64 | carry_bit;
    input[2] = (input1_u128 >> 1) as u64;
    input[3] = (input1_u128 >> 65) as u64;
    input
}
