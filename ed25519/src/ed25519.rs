use alloc::vec::Vec;
use core::fmt::{Debug, Display, Formatter};
use core::hash::{Hash, Hasher};
use core::iter::{Product, Sum};
use core::mem::transmute;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use core::{array, fmt};

use num_bigint::BigUint;
use p3_field::integers::QuotientMap;
use p3_field::op_assign_macros::{
    impl_add_assign, impl_div_methods, impl_mul_methods, impl_sub_assign, ring_sum,
};
use p3_field::{
    Field, InjectiveMonomial, Packable, PermutationMonomial, PrimeCharacteristicRing, PrimeField,
    RawDataSerializable, quotient_map_small_int,
};
use rand::distr::{Distribution, StandardUniform};
use rand::{Rng, RngExt};
use serde::{Deserialize, Deserializer, Serialize};

use crate::helpers::{halve_ed25519, monty_mul, to_biguint, wrapping_add, wrapping_sub};

/// The ed25519 base field prime `p = 2^255 - 19` as a little-endian array of 4-u64s.
///
/// Equal to: `57896044618658097711785492504343953926634992332820282019728792003956564819949`
pub(crate) const ED25519_PRIME: [u64; 4] = [
    0xffffffffffffffed,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x7fffffffffffffff,
];

// We use the Montgomery representation of the ed25519 prime, with respect to the
// constant `R = 2^256`. Note that `R = 2^256 mod p = 38` and `R^2 mod p = 1444`,
// which are unusually small (a consequence of `p` being very close to `2^255`).

/// The value `P^{-1} mod 2^64` where `P` is the ed25519 prime.
pub(crate) const ED25519_MONTY_MU_64: u64 = 0x79435e50d79435e5;

/// The square of the Montgomery constant `R = 2^256 mod P` for the ed25519 field.
///
/// Elements of the field are represented in Montgomery form, by `aR mod P`.
/// This constant is equal to `R^2 mod P` and is useful for converting elements into Montgomery form.
///
/// Equal to: `1444`
pub(crate) const ED25519_MONTY_R_SQ: [u64; 4] = [0x00000000000005a4, 0, 0, 0];

/// The ed25519 base field, defined as `F_P` where `P = 2^255 - 19`.
///
/// **Alignment on zkvm.** OpenVM's modular precompile operates on
/// `#[repr(C, align(32))] struct Ed25519Fp([u8; 32])`, so an unaligned
/// `[u64; 4]` has to be *copied* into that type and back around every
/// operation — three 32-byte copies per multiply, measured at **178
/// instructions against 34 for the bare chip op**. Giving the struct the
/// chip's alignment makes the two layouts identical (both are 32 canonical
/// little-endian bytes), so [`zk::as_fp`] is a pointer cast and the copies
/// disappear. Costs nothing: the struct is already exactly 32 bytes, so
/// alignment adds no padding and arrays keep a stride of 32.
#[derive(Copy, Clone, Default, Eq, PartialEq)]
#[cfg_attr(target_os = "zkvm", repr(C, align(32)))]
#[must_use]
pub struct Ed25519 {
    /// The MONTY form of the field element, a 255-bit integer less than `P` saved as a collection of u64's using a little-endian order.
    pub(crate) value: [u64; 4],
}

impl Ed25519 {
    /// Create a new field element from any `[u64; 4]` in little-endian limb order.
    ///
    /// Any value is accepted and automatically reduced modulo P.
    #[inline]
    pub const fn new(value: [u64; 4]) -> Self {
        #[cfg(not(target_os = "zkvm"))]
        { Self::new_monty(monty_mul(ED25519_MONTY_R_SQ, value)) }
        #[cfg(target_os = "zkvm")]
        {
            // canonical reduce: value < 2^256 < 3P, so at most two subtractions of P.
            let (v1, u1) = wrapping_sub(value, ED25519_PRIME);
            let value = if u1 { value } else { v1 };
            let (v2, u2) = wrapping_sub(value, ED25519_PRIME);
            let value = if u2 { value } else { v2 };
            Self::new_monty(value)
        }
    }

    /// Convert a `[[u64; 4]; N]` array to an array of field elements.
    ///
    /// Const version of `input.map(Ed25519::new)`.
    #[inline]
    pub const fn new_array<const N: usize>(input: [[u64; 4]; N]) -> [Self; N] {
        let mut output = [Self::ZERO; N];
        let mut i = 0;
        while i < N {
            output[i] = Self::new(input[i]);
            i += 1;
        }
        output
    }

    /// Creates a new ed25519 field element from an array of 4 u64's.
    ///
    /// The array is assumed to correspond to a 255-bit integer less than P and is interpreted as
    /// already being in Montgomery form.
    #[inline]
    pub(crate) const fn new_monty(value: [u64; 4]) -> Self {
        Self { value }
    }

    #[inline]
    #[allow(clippy::needless_pass_by_value)]
    /// Construct from the canonical little-endian 4-word representation.
    ///
    /// The inverse of [`PrimeField::canonical_u64_words`]. On a zkVM guest the
    /// canonical form *is* the internal representation, so this is a copy —
    /// which is what makes a table of precomputed field elements readable for
    /// ~34 instructions instead of ~6,809 to re-derive one.
    ///
    /// # Panics
    ///
    /// Never, but values `>= P` are reduced on the host and taken as-is on the
    /// guest; pass canonical words.
    pub const fn from_canonical_words(words: [u64; 4]) -> Self {
        #[cfg(not(target_os = "zkvm"))]
        {
            Self::new_monty(monty_mul(ED25519_MONTY_R_SQ, words))
        }
        #[cfg(target_os = "zkvm")]
        {
            Self::new_monty(words)
        }
    }

    pub fn from_biguint(value: BigUint) -> Option<Self> {
        #[cfg(target_os = "zkvm")]
        { Some(Self::new_monty(zk::from_bu(value))) }
        #[cfg(not(target_os = "zkvm"))]
        {
            let digits = value.to_u64_digits();
            let num_dig = digits.len();
            match num_dig {
                0 => Some(Self::ZERO),
                1..=4 => {
                    let mut inner = [0; 4];
                    inner[..num_dig].copy_from_slice(&digits);
                    // A 4-limb value can be up to `2^256 - 1 > P`; `monty_mul(R^2, inner)` yields
                    // the Montgomery form of `inner mod P`, correctly reduced.
                    Some(Self::new_monty(monty_mul(ED25519_MONTY_R_SQ, inner)))
                }
                _ => None, // Too many digits for ed25519
            }
        }
    }

    /// Converts a byte array in little-endian order to a field element.
    ///
    /// Assumes the bytes correspond to the Montgomery form of the desired field element.
    ///
    /// Returns None if the byte array is not exactly 32 bytes long or if the value
    /// represented by the byte array is not less than the ed25519 prime.
    ///
    /// The monty-byte parser (inverse of [`Self::into_bytes`]); kept as a utility even though
    /// serde now uses the canonical little-endian encoding (portable host↔zkvm).
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn from_bytes_monty(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 32 {
            return None;
        }
        let value: [u64; 4] = array::from_fn(|i| {
            // Convert each 8 bytes to a u64 in little-endian order.
            let start = i * 8;
            let end = start + 8;
            // This unwrap is safe due to the length check above.
            u64::from_le_bytes(bytes[start..end].try_into().unwrap())
        });
        // Check if the value is less than the prime.
        if value.iter().rev().cmp(ED25519_PRIME.iter().rev()) == core::cmp::Ordering::Less {
            Some(Self::new_monty(value))
        } else {
            None
        }
    }
}

#[cfg(target_os = "zkvm")]
pub(crate) mod zk {
    use num_bigint::BigUint;
    use openvm_algebra_guest::IntMod;
    openvm_algebra_moduli_macros::moduli_declare! {
        Ed25519Fp { modulus = "57896044618658097711785492504343953926634992332820282019728792003956564819949" },
    }
    #[inline]
    /// Borrow the canonical words *as* the precompile's type, with no copy.
    ///
    /// # Safety
    ///
    /// `Ed25519Fp` is `#[repr(C, align(32))]` over `[u8; 32]` and the caller's
    /// `[u64; 4]` lives inside an `Ed25519`, which carries `align(32)` on this
    /// target — so size and alignment match. The byte patterns agree because
    /// both are the canonical value little-endian, which on zkvm is exactly
    /// what `Ed25519::value` holds.
    #[inline]
    pub(crate) fn as_fp(v: &[u64; 4]) -> &Ed25519Fp {
        const {
            assert!(core::mem::size_of::<Ed25519Fp>() == 32);
            assert!(core::mem::align_of::<Ed25519Fp>() == 32);
        }
        unsafe { &*core::ptr::from_ref(v).cast::<Ed25519Fp>() }
    }

    /// Inverse of [`as_fp`], by value — a move, so alignment is irrelevant.
    #[inline]
    pub(crate) fn into_words(x: Ed25519Fp) -> [u64; 4] {
        const {
            assert!(core::mem::size_of::<Ed25519Fp>() == 32);
        }
        unsafe { core::mem::transmute::<Ed25519Fp, [u64; 4]>(x) }
    }

    #[allow(dead_code)]
    fn to_fp(v: [u64; 4]) -> Ed25519Fp {
        let bytes: [u8; 32] = unsafe { core::mem::transmute(v.map(u64::to_le_bytes)) };
        Ed25519Fp::from_le_bytes_unchecked(&bytes)
    }
    #[inline]
    #[allow(dead_code)]
    fn from_fp(x: &Ed25519Fp) -> [u64; 4] {
        let b = x.as_le_bytes();
        core::array::from_fn(|i| u64::from_le_bytes(b[i * 8..i * 8 + 8].try_into().unwrap()))
    }
    #[inline]
    pub(crate) fn mul(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
        into_words(as_fp(a) * as_fp(b))
    }
    #[inline]
    pub(crate) fn from_bu(b: BigUint) -> [u64; 4] {
        // b is always < p (curve constants + validated deserialize); canonical repr is
        // just its little-endian bytes packed into [u64; 4] (no IntMod/num-bigint-version dep).
        let bytes = b.to_bytes_le();
        let mut buf = [0u8; 32];
        let n = core::cmp::min(bytes.len(), 32);
        buf[..n].copy_from_slice(&bytes[..n]);
        core::array::from_fn(|i| u64::from_le_bytes(buf[i * 8..i * 8 + 8].try_into().unwrap()))
    }
}

impl Ed25519 {
    /// Parse the canonical little-endian 32-byte form, rejecting anything
    /// `>= P`.
    ///
    /// No `BigUint`: the bytes become four `u64`s, the bound check is a
    /// four-word comparison, and construction is [`Self::from_canonical_words`].
    /// The obvious implementation allocated a `Vec<u8>`, two `BigUint`s **and a
    /// third for the modulus on every call**, which cost ~1,830 instructions
    /// per element in a zkVM guest — 57% of the Quasar recursion guest was
    /// deserialising its own input.
    #[must_use]
    pub fn from_canonical_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 32 {
            return None;
        }
        let words: [u64; 4] = core::array::from_fn(|i| {
            u64::from_le_bytes(
                bytes[i * 8..i * 8 + 8]
                    .try_into()
                    .expect("eight bytes"),
            )
        });
        // words < P, little-endian: compare from the most significant word down.
        let mut canonical = false;
        for i in (0..4).rev() {
            if words[i] < ED25519_PRIME[i] {
                canonical = true;
                break;
            }
            if words[i] > ED25519_PRIME[i] {
                return None;
            }
        }
        if !canonical {
            return None; // equal to P
        }
        Some(Self::from_canonical_words(words))
    }

    /// The canonical little-endian 32-byte form.
    #[must_use]
    pub fn to_canonical_bytes(&self) -> [u8; 32] {
        let mut words = [0u64; 8];
        let _ = <Self as p3_field::PrimeField>::canonical_u64_words(self, &mut words);
        let mut bytes = [0u8; 32];
        for i in 0..4 {
            bytes[i * 8..i * 8 + 8].copy_from_slice(&words[i].to_le_bytes());
        }
        bytes
    }
}

impl Serialize for Ed25519 {
    /// Serializes to the canonical little-endian 32-byte representation (matches host & zkvm reprs).
    #[inline]
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_bytes(&self.to_canonical_bytes())
    }
}

impl<'de> Deserialize<'de> for Ed25519 {
    /// Deserializes the canonical little-endian 32-byte form, rejecting values
    /// `>= P` (so an encoding is unique and the transcript unambiguous).
    #[inline]
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct Visitor;

        impl<'de> serde::de::Visitor<'de> for Visitor {
            type Value = Ed25519;

            fn expecting(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                f.write_str("32 canonical little-endian bytes")
            }

            fn visit_bytes<E: serde::de::Error>(self, bytes: &[u8]) -> Result<Ed25519, E> {
                Ed25519::from_canonical_bytes(bytes)
                    .ok_or_else(|| E::custom("Invalid field element"))
            }

            fn visit_borrowed_bytes<E: serde::de::Error>(
                self,
                bytes: &'de [u8],
            ) -> Result<Ed25519, E> {
                self.visit_bytes(bytes)
            }

            /// Formats that present bytes as a sequence (JSON, say).
            ///
            /// A wrong length is reported as `Invalid field element` rather than
            /// serde's `invalid_length`, so every rejection reason a caller can
            /// trigger reads the same — the length check and the `>= P` check
            /// are both "this is not a field element".
            fn visit_seq<A: serde::de::SeqAccess<'de>>(
                self,
                mut seq: A,
            ) -> Result<Ed25519, A::Error> {
                use serde::de::Error;
                let invalid = || A::Error::custom("Invalid field element");
                let mut buf = [0u8; 32];
                for slot in &mut buf {
                    *slot = seq.next_element()?.ok_or_else(invalid)?;
                }
                // Too long is as wrong as too short.
                if seq.next_element::<u8>()?.is_some() {
                    return Err(invalid());
                }
                Ed25519::from_canonical_bytes(&buf).ok_or_else(invalid)
            }
        }

        d.deserialize_bytes(Visitor)
    }
}

impl Packable for Ed25519 {}

impl Hash for Ed25519 {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        for limb in self.value.as_ref() {
            state.write_u64(*limb);
        }
    }
}

impl Ord for Ed25519 {
    #[inline]
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.value.iter().rev().cmp(other.value.iter().rev())
    }
}

impl PartialOrd for Ed25519 {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for Ed25519 {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        core::fmt::Display::fmt(&self.as_canonical_biguint(), f)
    }
}

impl Debug for Ed25519 {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        core::fmt::Debug::fmt(&self.as_canonical_biguint(), f)
    }
}

impl PrimeCharacteristicRing for Ed25519 {
    type PrimeSubfield = Self;

    const ZERO: Self = Self::new_monty([0, 0, 0, 0]);

    /// The Montgomery form of the ed25519 field element 1.
    ///
    /// Equal to `2^256 mod P = 38`.
    #[cfg(not(target_os = "zkvm"))]
    const ONE: Self = Self::new_monty([0x0000000000000026, 0, 0, 0]);
    #[cfg(target_os = "zkvm")]
    const ONE: Self = Self::new_monty([1, 0, 0, 0]);

    /// The Montgomery form of the ed25519 field element 2.
    ///
    /// Equal to `2^257 mod P = 76`.
    #[cfg(not(target_os = "zkvm"))]
    const TWO: Self = Self::new_monty([0x000000000000004c, 0, 0, 0]);
    #[cfg(target_os = "zkvm")]
    const TWO: Self = Self::new_monty([2, 0, 0, 0]);

    /// The Montgomery form of the ed25519 field element -1.
    ///
    /// Equal to `-2^256 mod P = P - 38 = 57896044618658097711785492504343953926634992332820282019728792003956564819911`.
    #[cfg(not(target_os = "zkvm"))]
    const NEG_ONE: Self = Self::new_monty([
        0xffffffffffffffc7,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
    ]);
    #[cfg(target_os = "zkvm")]
    const NEG_ONE: Self = Self::new_monty([
        0xffffffffffffffec,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
    ]);

    #[inline]
    fn from_prime_subfield(f: Self::PrimeSubfield) -> Self {
        f
    }

    #[inline]
    fn halve(&self) -> Self {
        Self::new_monty(halve_ed25519(self.value))
    }

    /// Squaring, counted separately from `Mul` for profiling. Computes the product directly
    /// (does NOT route through `Mul`), so `square()` calls increment `SQ_TOTAL`, not `MUL_TOTAL`.
    #[inline]
    fn square(&self) -> Self {
        #[cfg(feature = "count-ops")] crate::counters::on_sq();
        #[cfg(not(target_os = "zkvm"))]
        { Self::new_monty(monty_mul(self.value, self.value)) }
        #[cfg(target_os = "zkvm")]
        { Self::new_monty(zk::mul(&self.value, &self.value)) }
    }
}

/// Degree of the smallest permutation polynomial for the ed25519 base field.
///
/// As `p - 1 = 2^2 * (2^253 - 5)` is divisible by 2 and by 3 (`p - 1 = 2^255 - 20 ≡ 0 mod 3`),
/// the smallest choice for a degree `D` satisfying `gcd(p - 1, D) = 1` is 5.
impl InjectiveMonomial<5> for Ed25519 {}

/// The value `5^{-1} mod (P - 1)` where `P` is the ed25519 prime.
///
/// Satisfies `5 * ED25519_FIFTH_ROOT_EXP = 1 mod (P - 1)`.
const ED25519_FIFTH_ROOT_EXP: [u64; 4] = [
    0xccccccccccccccc1,
    0xcccccccccccccccc,
    0xcccccccccccccccc,
    0x4ccccccccccccccc,
];

/// The exponent `P - 2` used for Fermat inversion (`a^{P-2} = a^{-1}`).
const ED25519_PRIME_MINUS_2: [u64; 4] = [
    0xffffffffffffffeb,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x7fffffffffffffff,
];

/// Compute `x -> x^exp` using 4-bit windowed square-and-multiply.
fn exp_u256(base: Ed25519, exp: [u64; 4]) -> Ed25519 {
    let mut table = [Ed25519::ONE; 16];
    table[1] = base;
    for i in 2..16 {
        table[i] = table[i - 1] * base;
    }

    let mut started = false;
    let mut result = Ed25519::ONE;

    for nibble_idx in (0..64).rev() {
        let limb_idx = nibble_idx / 16;
        let shift = (nibble_idx % 16) * 4;
        let nibble = ((exp[limb_idx] >> shift) & 0xF) as usize;

        if started {
            result = result.exp_power_of_2(4);
            if nibble != 0 {
                result *= table[nibble];
            }
        } else if nibble != 0 {
            result = table[nibble];
            started = true;
        }
    }
    result
}

impl PermutationMonomial<5> for Ed25519 {
    /// In the field `Ed25519`, `a^{1/5}` is equal to `a^{ED25519_FIFTH_ROOT_EXP}`.
    ///
    /// This follows from the calculation `5 * ED25519_FIFTH_ROOT_EXP = 1 mod (P - 1)`.
    fn injective_exp_root_n(&self) -> Self {
        exp_u256(*self, ED25519_FIFTH_ROOT_EXP)
    }
}

impl RawDataSerializable for Ed25519 {
    const NUM_BYTES: usize = 32;

    #[allow(refining_impl_trait)]
    #[inline]
    fn into_bytes(self) -> [u8; 32] {
        // The transmute here maps from [[u8; 8]; 4] to [u8; 32] so is clearly safe.
        unsafe { transmute(self.value.map(|x| x.to_le_bytes())) }
    }

    #[inline]
    fn into_u32_stream(input: impl IntoIterator<Item = Self>) -> impl IntoIterator<Item = u32> {
        input.into_iter().flat_map(|x| {
            x.value
                .into_iter()
                .flat_map(|digit| [digit as u32, (digit >> 32) as u32])
        })
    }

    #[inline]
    fn into_u64_stream(input: impl IntoIterator<Item = Self>) -> impl IntoIterator<Item = u64> {
        input.into_iter().flat_map(|x| x.value)
    }

    #[inline]
    fn into_parallel_byte_streams<const N: usize>(
        input: impl IntoIterator<Item = [Self; N]>,
    ) -> impl IntoIterator<Item = [u8; N]> {
        input.into_iter().flat_map(|vector| {
            let bytes = vector.map(|elem| elem.into_bytes());
            (0..Self::NUM_BYTES).map(move |i| array::from_fn(|j| bytes[j][i]))
        })
    }

    #[inline]
    fn into_parallel_u32_streams<const N: usize>(
        input: impl IntoIterator<Item = [Self; N]>,
    ) -> impl IntoIterator<Item = [u32; N]> {
        input.into_iter().flat_map(|vector| {
            // The transmute here maps from [[u32; 2]; 4] to [u32; 8] so is clearly safe.
            let u32s: [[u32; 8]; N] = vector
                .map(|elem| unsafe { transmute(elem.value.map(|x| [x as u32, (x >> 32) as u32])) });
            (0..(Self::NUM_BYTES / 4)).map(move |i| array::from_fn(|j| u32s[j][i]))
        })
    }

    #[inline]
    fn into_parallel_u64_streams<const N: usize>(
        input: impl IntoIterator<Item = [Self; N]>,
    ) -> impl IntoIterator<Item = [u64; N]> {
        input.into_iter().flat_map(|vector| {
            let u64s = vector.map(|elem| elem.value);
            (0..(Self::NUM_BYTES / 8)).map(move |i| array::from_fn(|j| u64s[j][i]))
        })
    }
}

impl Field for Ed25519 {
    type Packing = Self;

    /// A plausible multiplicative generator for the ed25519 field: the canonical value `2`.
    ///
    /// **WARNING: this is NOT verified to be a true multiplicative generator.**
    /// `2` is the smallest quadratic non-residue mod `P` (so it generates the full 2-power part
    /// of the group), but confirming it generates the entire group of order
    /// `P - 1 = 2^2 * (2^253 - 5)` would require the prime factorization of the 253-bit factor
    /// `2^253 - 5`, which is infeasible.
    ///
    /// This is acceptable here because the intended use (a Brakedown-based STARK using this field
    /// as both trace and challenge field) never uses roots of unity or the generator: this field
    /// is non-two-adic and does not implement `TwoAdicField`. `GENERATOR` exists only to satisfy
    /// the `Field` trait bound.
    ///
    /// Equal to the Montgomery form of `2`, i.e. `2 * 2^256 mod P = 76`.
    #[cfg(not(target_os = "zkvm"))]
    const GENERATOR: Self = Self::new_monty([0x000000000000004c, 0, 0, 0]);
    #[cfg(target_os = "zkvm")]
    const GENERATOR: Self = Self::new_monty([2, 0, 0, 0]);

    #[inline]
    fn is_zero(&self) -> bool {
        self.value.iter().all(|&x| x == 0)
    }

    #[inline]
    fn try_inverse(&self) -> Option<Self> {
        // Fermat inversion: for a != 0, `a^{-1} = a^{P-2} mod P`.
        //
        // We use exponentiation-based inversion (rather than porting bn254's bit-length-specific
        // binary-GCD inversion, whose round counts and adjustment constants are tuned for the
        // 254-bit bn254 prime) as it is unconditionally correct for this 255-bit prime and reuses
        // the existing windowed exponentiation. Inversion is rare and can be batched.
        #[cfg(feature = "count-ops")] crate::counters::on_inv();
        (!self.is_zero()).then(|| exp_u256(*self, ED25519_PRIME_MINUS_2))
    }

    /// `r = 57896044618658097711785492504343953926634992332820282019728792003956564819949`
    #[inline]
    fn order() -> BigUint {
        to_biguint(ED25519_PRIME)
    }
}

quotient_map_small_int!(Ed25519, u128, [u8, u16, u32, u64]);
quotient_map_small_int!(Ed25519, i128, [i8, i16, i32, i64]);

impl QuotientMap<u128> for Ed25519 {
    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    fn from_int(int: u128) -> Self {
        // Need to convert into Monty form. As the monty reduction strips out a factor of `R`,
        // we can do this by multiplying by `R^2` and doing a monty reduction.
        #[cfg(not(target_os = "zkvm"))]
        { let monty_form = monty_mul(ED25519_MONTY_R_SQ, [int as u64, (int >> 64) as u64, 0, 0]); Self::new_monty(monty_form) }
        #[cfg(target_os = "zkvm")]
        { Self::new_monty([int as u64, (int >> 64) as u64, 0, 0]) }
    }

    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    fn from_canonical_checked(int: u128) -> Option<Self> {
        Some(Self::from_int(int))
    }

    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    unsafe fn from_canonical_unchecked(int: u128) -> Self {
        Self::from_int(int)
    }
}

impl QuotientMap<i128> for Ed25519 {
    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    fn from_int(int: i128) -> Self {
        // Nothing better than just branching based on the sign of int.
        if int >= 0 {
            Self::from_int(int as u128)
        } else {
            -Self::from_int((-int) as u128)
        }
    }

    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    fn from_canonical_checked(int: i128) -> Option<Self> {
        Some(Self::from_int(int))
    }

    /// Due to the size of the ed25519 prime, the input value is always canonical.
    #[inline]
    unsafe fn from_canonical_unchecked(int: i128) -> Self {
        Self::from_int(int)
    }
}

impl PrimeField for Ed25519 {
    #[inline]
    fn canonical_u64_words(&self, out: &mut [u64]) -> usize {
        // ed25519 is a fixed 4-word field, so the canonical value is already a
        // `[u64; 4]` — no `BigUint`, and in particular none of the byte-vector
        // round trip `to_biguint` does. This runs once per hashed element in the
        // cross-field Merkle/transcript path, which is the hottest loop in the
        // OpenVM recursion verifier.
        #[cfg(not(target_os = "zkvm"))]
        let value = monty_mul(self.value, [1, 0, 0, 0]);
        #[cfg(target_os = "zkvm")]
        let value = self.value;
        out[..4].copy_from_slice(&value);
        4
    }

    fn as_canonical_biguint(&self) -> BigUint {
        // `monty_mul` strips out a factor of `R` so multiplying by `1` converts a montgomery
        // representation into a canonical representation.
        #[cfg(not(target_os = "zkvm"))]
        { let out_val = monty_mul(self.value, [1, 0, 0, 0]); to_biguint(out_val) }
        #[cfg(target_os = "zkvm")]
        { to_biguint(self.value) }
    }
}

impl Add for Ed25519 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        #[cfg(feature = "count-ops")] crate::counters::on_add();
        let lhs = self.value;
        let rhs = rhs.value;

        let (sum, overflow) = wrapping_add(lhs, rhs);

        // As the inputs are < 2^255, the output is < 2^256 so the final carry flag should be false.
        debug_assert!(!overflow);

        // If output is bigger than ED25519_PRIME, we should subtract ED25519_PRIME from it.
        let (sum_corr, underflow) = wrapping_sub(sum, ED25519_PRIME);

        if underflow {
            Self::new_monty(sum)
        } else {
            Self::new_monty(sum_corr)
        }
    }
}

impl Sub for Ed25519 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        #[cfg(feature = "count-ops")] crate::counters::on_add();
        let lhs = self.value;
        let rhs = rhs.value;

        let (mut sub, underflow) = wrapping_sub(lhs, rhs);

        // if lhs < rhs, we need to add ED25519_PRIME to the result.
        if underflow {
            (sub, _) = wrapping_add(sub, ED25519_PRIME);
        }

        Self::new_monty(sub)
    }
}

impl Neg for Ed25519 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::ZERO - self
    }
}

impl Mul for Ed25519 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        #[cfg(feature = "count-ops")] crate::counters::on_mul();
        #[cfg(not(target_os = "zkvm"))]
        { Self::new_monty(monty_mul(self.value, rhs.value)) }
        #[cfg(target_os = "zkvm")]
        { Self::new_monty(zk::mul(&self.value, &rhs.value)) }
    }
}

impl_add_assign!(Ed25519);
impl_sub_assign!(Ed25519);
impl_mul_methods!(Ed25519);
ring_sum!(Ed25519);
impl_div_methods!(Ed25519, Ed25519);

impl Distribution<Ed25519> for StandardUniform {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Ed25519 {
        // Rejection sampling constructed BY VALUE (repr-independent) so the same RNG stream yields
        // the same field element on the Montgomery host and the canonical zkvm guest
        // (native-field-fork). `from_bytes_monty` treats the bytes as internal limbs (Montgomery vs
        // canonical) -> different values, which breaks cross-platform Brakedown/Basefold code matrices.
        // Value-identical to `from_biguint` on a masked 32-byte draw (which has
        // at most four u64 digits, so that path never returns `None`), but with
        // no `BigUint`: the round trip through `from_bytes_le` and `to_bytes_le`
        // allocated twice per element and measured **5,127 instructions in an
        // rv32 zkVM guest**, against ~190 for a field multiplication. Building a
        // foldable code samples one of these per twist, so it was the dominant
        // term in the OpenVM recursion guest's setup.
        let mut trial_element: [u8; 32] = rng.random();
        trial_element[31] &= (1_u8 << 7) - 1;
        let inner: [u64; 4] = core::array::from_fn(|i| {
            u64::from_le_bytes(
                trial_element[i * 8..i * 8 + 8]
                    .try_into()
                    .expect("eight bytes"),
            )
        });
        // A 4-limb value can exceed P; `monty_mul(R², inner)` is `inner mod P`
        // in Montgomery form. The guest's representation is already canonical.
        #[cfg(not(target_os = "zkvm"))]
        {
            Ed25519::new_monty(monty_mul(ED25519_MONTY_R_SQ, inner))
        }
        #[cfg(target_os = "zkvm")]
        {
            Ed25519::new_monty(inner)
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use p3_field_testing::{test_field_json_serialization, test_prime_field};
    use proptest::prelude::*;

    use super::*;

    type F = Ed25519;

    #[test]
    fn test_ed25519() {
        let big_int_100 = BigUint::from(100u32);
        let big_int_p = to_biguint(ED25519_PRIME);
        let big_int_2_256_min_1 = to_biguint([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
        ]);
        // `2^256 mod P = 38`.
        let big_int_2_256_mod_p = BigUint::from(38u32);

        let f_100 = F::from_biguint(big_int_100.clone()).unwrap();
        assert_eq!(f_100.as_canonical_biguint(), BigUint::from(100u32));
        assert_eq!(F::from_biguint(BigUint::ZERO), Some(F::ZERO));

        // Unlike bn254 (whose prime is ~253.6 bits, so `5P < 2^256`), here `P ~ 2^255`, so only
        // multiples `k*P` with `k*P < 2^256` (i.e. `k <= 2`) fit into `from_biguint`'s 4-limb input.
        // Larger multiples correctly return `None` (more than 4 u64 digits), tested separately below.
        for i in 0_u32..3_u32 {
            assert_eq!(F::from_biguint(big_int_p.clone() * i), Some(F::ZERO));
        }
        // `100 + P < 2^256` fits and must reduce to 100.
        assert_eq!(
            F::from_biguint(big_int_100 + big_int_p.clone()),
            Some(f_100)
        );
        // `3*P > 2^256` and `6*P > 2^256`: too many digits, rejected.
        assert_eq!(F::from_biguint(big_int_p.clone() * 3_u32), None);
        assert_eq!(F::from_biguint(big_int_p * 6_u32), None);
        assert_eq!(
            F::from_biguint(big_int_2_256_min_1).unwrap(),
            F::NEG_ONE + F::from_biguint(big_int_2_256_mod_p).unwrap()
        );

        // Basic Montgomery constant sanity checks.
        assert_eq!(F::ZERO.as_canonical_biguint(), BigUint::ZERO);
        assert_eq!(F::ONE.as_canonical_biguint(), BigUint::from(1u32));
        assert_eq!(F::TWO.as_canonical_biguint(), BigUint::from(2u32));
        assert_eq!(F::ONE + F::NEG_ONE, F::ZERO);
        assert_eq!(F::NEG_ONE.as_canonical_biguint(), to_biguint(ED25519_PRIME) - 1u32);

        // Documented (unverified) generator has canonical value 2.
        assert_eq!(F::GENERATOR, F::from_u8(2));
        assert_eq!(F::GENERATOR.as_canonical_biguint(), BigUint::from(2u32));

        let f_1 = F::ONE;
        let f_2 = F::TWO;
        let f_r_minus_1 = F::NEG_ONE;
        let f_r_minus_2 = F::NEG_ONE + F::NEG_ONE;

        test_field_json_serialization(&[f_100, f_1, f_2, f_r_minus_1, f_r_minus_2]);
    }

    #[test]
    fn test_add_sub_neg_roundtrips() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;
        let mut rng = SmallRng::seed_from_u64(42);
        for _ in 0..1000 {
            let a: F = rng.random();
            let b: F = rng.random();
            // (a + b) - b == a
            assert_eq!((a + b) - b, a);
            // a + (-a) == 0
            assert_eq!(a + (-a), F::ZERO);
            // -(-a) == a
            assert_eq!(-(-a), a);
            // a - b == a + (-b)
            assert_eq!(a - b, a + (-b));
            // a.halve().double() == a
            assert_eq!(a.halve().double(), a);
        }
    }

    #[test]
    fn test_mul_inverse_roundtrip() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;
        let mut rng = SmallRng::seed_from_u64(7);
        for _ in 0..500 {
            let a: F = rng.random();
            if a.is_zero() {
                continue;
            }
            let inv = a.try_inverse().unwrap();
            assert_eq!(a * inv, F::ONE);
            assert_eq!(a.inverse(), inv);
        }
        assert_eq!(F::ZERO.try_inverse(), None);
    }

    #[test]
    fn test_from_int_canonical_roundtrip() {
        for v in [0u128, 1, 2, 5, 38, 76, 12345, u64::MAX as u128, u128::MAX] {
            let f = F::from_int(v);
            let expected = BigUint::from(v) % to_biguint(ED25519_PRIME);
            assert_eq!(f.as_canonical_biguint(), expected);
        }
    }

    #[test]
    fn test_deserialize_boundary_rejections() {
        // Valid inputs are exactly the 32-byte little-endian encodings of values in [0, P).

        // value == P.
        let prime_bytes: [u8; 32] = unsafe { transmute(ED25519_PRIME.map(u64::to_le_bytes)) };

        // value == P + 1. Low limb is 0xffffffffffffffed, so `+= 1` cannot overflow.
        let mut over_prime = ED25519_PRIME;
        over_prime[0] += 1;
        let over_prime_bytes: [u8; 32] = unsafe { transmute(over_prime.map(u64::to_le_bytes)) };

        let cases: &[&[u8]] = &[
            &prime_bytes,      // == P
            &over_prime_bytes, // == P + 1
            &[0xff_u8; 32],    // == 2^256 - 1, top of the byte-space
            &[],               // empty: fails the length check before any limb is read
            &[0u8; 31],        // one byte too short
            &[0u8; 33],        // one byte too long
        ];

        for bytes in cases {
            let json = serde_json::to_string(bytes).unwrap();
            let err = serde_json::from_str::<Ed25519>(&json).unwrap_err();
            assert!(
                err.to_string().contains("Invalid field element"),
                "expected field-validation error, got: {err}"
            );
        }
    }

    #[test]
    fn test_deserialize_accepts_max_valid() {
        // The boundary at P is exclusive: P - 1 must still round-trip.
        // PRIME's low limb is non-zero, so `-= 1` cannot underflow.
        let mut p_minus_one = ED25519_PRIME;
        p_minus_one[0] -= 1;
        let bytes: [u8; 32] = unsafe { transmute(p_minus_one.map(u64::to_le_bytes)) };

        let json = serde_json::to_string(bytes.as_slice()).unwrap();
        let _: Ed25519 = serde_json::from_str(&json).expect("P-1 must round-trip");
    }

    // `test_prime_field!` only exercises small-integer `QuotientMap` conversions; it does NOT
    // require the prime factorization of `P - 1` (that is only needed by `test_field!`'s generator
    // check, which we intentionally omit since factoring `2^253 - 5` is infeasible).
    test_prime_field!(crate::Ed25519);

    fn arb_ed25519() -> impl Strategy<Value = F> {
        prop::array::uniform4(prop::num::u64::ANY).prop_map(F::new)
    }

    proptest! {
        // Validates the N=5 permutation and the `ED25519_FIFTH_ROOT_EXP = 5^{-1} mod (P-1)`
        // constant, exercising the full field arithmetic (Montgomery mul + exponentiation).
        #[test]
        fn test_permutation_monomial_5(x in arb_ed25519()) {
            assert_eq!(x.injective_exp_n().injective_exp_root_n(), x);
            assert_eq!(x.injective_exp_root_n().injective_exp_n(), x);
        }
    }
}
