//! A Poseidon2 permutation over the ed25519 base field (`p = 2^255 - 19`).
//!
//! This mirrors the bn254 Poseidon2 wiring (`p3_bn254::poseidon2`): a width-3
//! state with S-box degree 5 (the field implements [`InjectiveMonomial<5>`]),
//! the `1 + Diag([1, 1, 2])` internal diffusion matrix, and the standard
//! external MDS layer.
//!
//! # DEMO-GRADE CONSTANTS
//!
//! Unlike bn254 (which imports audited round constants from the HorizenLabs
//! reference), there is **no** published Poseidon2 instance for the ed25519 base
//! field. The round constants here are therefore generated from a seeded RNG via
//! [`Poseidon2::new_from_rng`]. This is sufficient for an end-to-end
//! demonstration and testing of the field-agnostic Brakedown STARK, but it is
//! **NOT production-secure**: the constants have not undergone any cryptographic
//! analysis. Do not use this permutation to secure real value.

use alloc::vec::Vec;

use p3_field::PrimeCharacteristicRing;
use p3_poseidon2::{
    ExternalLayer, ExternalLayerConstants, ExternalLayerConstructor, HLMDSMat4, InternalLayer,
    InternalLayerConstructor, Poseidon2, add_rc_and_sbox_generic, external_initial_permute_state,
    external_terminal_permute_state, internal_permute_state,
};

use crate::Ed25519;

/// S-box degree for ed25519 Poseidon2.
///
/// As `p - 1 = 2^2 * (2^253 - 5)` is divisible by 2 and 3, the smallest valid
/// exponent satisfying `gcd(α, p - 1) = 1` is 5 (matches [`InjectiveMonomial<5>`]).
pub const ED25519_S_BOX_DEGREE: u64 = 5;

/// Total number of full (external) rounds for ed25519 Poseidon2 at width 3.
///
/// Mirrors the bn254 width-3 instance (`R_F = 8`, i.e. 4 per side).
pub const ED25519_POSEIDON2_FULL_ROUNDS: usize = 8;

/// Number of partial (internal) rounds for ed25519 Poseidon2 at width 3.
///
/// Mirrors the bn254 width-3 instance (`R_P = 56`).
pub const ED25519_POSEIDON2_PARTIAL_ROUNDS_3: usize = 56;

/// The only supported state width for ed25519 Poseidon2 (matches bn254).
const ED25519_WIDTH: usize = 3;

/// An implementation of the Poseidon2 permutation for the [`Ed25519`] field.
///
/// It acts on arrays of the form `[Ed25519; WIDTH]`.
pub type Poseidon2Ed25519<const WIDTH: usize> = Poseidon2<
    Ed25519,
    Poseidon2ExternalLayerEd25519<WIDTH>,
    Poseidon2InternalLayerEd25519,
    WIDTH,
    ED25519_S_BOX_DEGREE,
>;

#[derive(Debug, Clone, Default)]
pub struct Poseidon2InternalLayerEd25519 {
    internal_constants: Vec<Ed25519>,
}

impl InternalLayerConstructor<Ed25519> for Poseidon2InternalLayerEd25519 {
    fn new_from_constants(internal_constants: Vec<Ed25519>) -> Self {
        Self { internal_constants }
    }
}

/// A faster version of `matmul_internal` making use of the fact that
/// the internal matrix is equal to:
/// ```ignore
///                             [2, 1, 1]
///     1 + Diag([1, 1, 2]) =   [1, 2, 1]
///                             [1, 1, 3]
/// ```
fn ed25519_matmul_internal(state: &mut [Ed25519; 3]) {
    let sum = state[0] + (state[1] + state[2]);
    state[0] += sum;
    state[1] += sum;
    state[2] = state[2].double() + sum;
}

impl InternalLayer<Ed25519, ED25519_WIDTH, ED25519_S_BOX_DEGREE> for Poseidon2InternalLayerEd25519 {
    fn permute_state(&self, state: &mut [Ed25519; ED25519_WIDTH]) {
        internal_permute_state(state, ed25519_matmul_internal, &self.internal_constants);
    }
}

pub type Poseidon2ExternalLayerEd25519<const WIDTH: usize> = ExternalLayerConstants<Ed25519, WIDTH>;

impl<const WIDTH: usize> ExternalLayerConstructor<Ed25519, WIDTH>
    for Poseidon2ExternalLayerEd25519<WIDTH>
{
    fn new_from_constants(external_constants: Self) -> Self {
        external_constants
    }
}

impl<const WIDTH: usize> ExternalLayer<Ed25519, WIDTH, ED25519_S_BOX_DEGREE>
    for Poseidon2ExternalLayerEd25519<WIDTH>
{
    fn permute_state_initial(&self, state: &mut [Ed25519; WIDTH]) {
        external_initial_permute_state(
            state,
            self.get_initial_constants(),
            add_rc_and_sbox_generic,
            &HLMDSMat4,
        );
    }

    fn permute_state_terminal(&self, state: &mut [Ed25519; WIDTH]) {
        external_terminal_permute_state(
            state,
            self.get_terminal_constants(),
            add_rc_and_sbox_generic,
            &HLMDSMat4,
        );
    }
}

#[cfg(test)]
mod tests {
    use p3_symmetric::Permutation;
    use rand::rngs::SmallRng;
    use rand::{RngExt, SeedableRng};

    use super::*;

    const WIDTH: usize = 3;

    fn build_perm(seed: u64) -> Poseidon2Ed25519<WIDTH> {
        let mut rng = SmallRng::seed_from_u64(seed);
        Poseidon2Ed25519::<WIDTH>::new_from_rng(
            ED25519_POSEIDON2_FULL_ROUNDS,
            ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
            &mut rng,
        )
    }

    #[test]
    fn poseidon2_runs_and_is_deterministic() {
        let perm = build_perm(0xED25519);
        let mut rng = SmallRng::seed_from_u64(1);
        let input = rng.random::<[Ed25519; WIDTH]>();

        let mut out_a = input;
        perm.permute_mut(&mut out_a);
        let mut out_b = input;
        perm.permute_mut(&mut out_b);

        // Deterministic: same permutation + same input => same output.
        assert_eq!(out_a, out_b);
        // The permutation actually mixes the state (not the identity).
        assert_ne!(out_a, input);
    }

    #[test]
    fn poseidon2_same_seed_same_permutation() {
        // Two permutations built from the same seed must agree on all inputs.
        let perm1 = build_perm(42);
        let perm2 = build_perm(42);
        let mut rng = SmallRng::seed_from_u64(7);
        for _ in 0..8 {
            let input = rng.random::<[Ed25519; WIDTH]>();
            let mut a = input;
            let mut b = input;
            perm1.permute_mut(&mut a);
            perm2.permute_mut(&mut b);
            assert_eq!(a, b);
        }
    }

    #[test]
    fn poseidon2_is_a_permutation() {
        // Distinct inputs map to distinct outputs (bijectivity spot-check).
        let perm = build_perm(99);
        let mut rng = SmallRng::seed_from_u64(3);
        let x = rng.random::<[Ed25519; WIDTH]>();
        let mut y = x;
        y[0] += Ed25519::ONE;

        let mut px = x;
        let mut py = y;
        perm.permute_mut(&mut px);
        perm.permute_mut(&mut py);
        assert_ne!(px, py);
    }
}
