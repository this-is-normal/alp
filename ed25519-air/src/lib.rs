//! AIRs for elliptic-curve arithmetic on the prime-order subgroup of ed25519.
//!
//! These AIRs assume the STARK's native field *is* the curve's base field (`2^255 - 19`),
//! i.e. no foreign-field / non-native-modulus emulation is performed here. That base-field
//! implementation is being developed elsewhere; until it lands, tests in this crate run over
//! an arbitrary native field (e.g. `p3-bn254`) supplied with placeholder curve constants, to
//! validate the constraint structure independently of the concrete field.
//!
//! Points are represented in affine `(x, y)` coordinates on the twisted Edwards curve
//! `a*x^2 + y^2 = 1 + d*x^2*y^2` with `a = -1`, using the complete addition law that is valid
//! for both distinct-point addition and doubling (no exceptional cases).

#![no_std]

extern crate alloc;

pub mod add;
pub mod eddsa;
pub mod fixed_scalarmul;
pub mod pointcheck;
pub mod scalarmul;
pub mod subgroupcheck;
