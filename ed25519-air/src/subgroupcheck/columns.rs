use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

use crate::add::NUM_EC_ADD_COLS;
use crate::pointcheck::NUM_EC_POINTCHECK_COLS;

/// Trace columns for a subgroup-membership certificate: prove that `(q_x, q_y)` is `8*P` for
/// some on-curve witness point `P`, without revealing `P` beyond what's needed to check it.
///
/// Since `E(F)` is cyclic of order `8*l` for ed25519 (cofactor 8, prime `l`), the image of the
/// multiplication-by-8 map is exactly the unique order-`l` subgroup — so exhibiting such a `P`
/// is equivalent to proving `(q_x, q_y)` lies in that prime-order subgroup. See
/// `docs/plans/SubgroupCheck.md` for the full argument.
///
/// Composed of four sub-air blocks tied together with cheap degree-1 equality constraints (see
/// [`crate::subgroupcheck::air::EcSubgroupCheckAir`]):
/// - `cc_*` — an [`crate::pointcheck::EcPointCheckAir`] sub-air proving `P = (p_x, p_y)` is
///   on-curve.
/// - `d1_*` — an [`crate::add::EcAddAir`] sub-air computing `2*P`.
/// - `d2_*` — computing `4*P` from `2*P`.
/// - `d3_*` — computing `8*P` from `4*P`; its output is constrained equal to `(q_x, q_y)`.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcSubgroupCheckCols<T> {
    /// Witness point `P` (private): any on-curve point with `8*P = Q`.
    pub p_x: T,
    pub p_y: T,

    // -- point-check sub-air block: x,y = copies of (p_x,p_y); x2,y2 = helpers --
    pub cc_x: T,
    pub cc_y: T,
    pub cc_x2: T,
    pub cc_y2: T,

    // -- doubling sub-air block 1: P + P = 2P --
    pub d1_x1: T,
    pub d1_y1: T,
    pub d1_x2: T,
    pub d1_y2: T,
    pub d1_x3: T,
    pub d1_y3: T,
    pub d1_p: T,
    pub d1_q: T,
    pub d1_m: T,

    // -- doubling sub-air block 2: 2P + 2P = 4P --
    pub d2_x1: T,
    pub d2_y1: T,
    pub d2_x2: T,
    pub d2_y2: T,
    pub d2_x3: T,
    pub d2_y3: T,
    pub d2_p: T,
    pub d2_q: T,
    pub d2_m: T,

    // -- doubling sub-air block 3: 4P + 4P = 8P --
    pub d3_x1: T,
    pub d3_y1: T,
    pub d3_x2: T,
    pub d3_y2: T,
    pub d3_x3: T,
    pub d3_y3: T,
    pub d3_p: T,
    pub d3_q: T,
    pub d3_m: T,

    /// The point under test: constrained equal to `8*P` (`d3_x3, d3_y3`).
    pub q_x: T,
    pub q_y: T,
}

pub const NUM_EC_SUBGROUPCHECK_COLS: usize = size_of::<EcSubgroupCheckCols<u8>>();

// The point-check and doubling sub-air blocks above must exactly match their AIRs' layout/width.
const _: () = assert!(4 == NUM_EC_POINTCHECK_COLS);
const _: () = assert!(9 == NUM_EC_ADD_COLS);

impl<T> Borrow<EcSubgroupCheckCols<T>> for [T] {
    fn borrow(&self) -> &EcSubgroupCheckCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_SUBGROUPCHECK_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EcSubgroupCheckCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcSubgroupCheckCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcSubgroupCheckCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_SUBGROUPCHECK_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EcSubgroupCheckCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

/// Column index (within the flat row) where the point-check sub-air block starts.
pub const CC_BLOCK_START: usize = 2;
/// Column index where the first doubling sub-air block (`P + P = 2P`) starts.
pub const D1_BLOCK_START: usize = CC_BLOCK_START + NUM_EC_POINTCHECK_COLS;
/// Column index where the second doubling sub-air block (`2P + 2P = 4P`) starts.
pub const D2_BLOCK_START: usize = D1_BLOCK_START + NUM_EC_ADD_COLS;
/// Column index where the third doubling sub-air block (`4P + 4P = 8P`) starts.
pub const D3_BLOCK_START: usize = D2_BLOCK_START + NUM_EC_ADD_COLS;
