use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Trace columns for checking that `(x, y)` lies on the twisted-Edwards curve
/// `a*x^2 + y^2 = 1 + d*x^2*y^2`.
///
/// `x2`, `y2` are helper witness columns that keep every constraint at degree <= 2:
/// - `x2 = x*x`
/// - `y2 = y*y`
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcPointCheckCols<T> {
    pub x: T,
    pub y: T,
    pub x2: T,
    pub y2: T,
}

pub const NUM_EC_POINTCHECK_COLS: usize = size_of::<EcPointCheckCols<u8>>();

impl<T> Borrow<EcPointCheckCols<T>> for [T] {
    fn borrow(&self) -> &EcPointCheckCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_POINTCHECK_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EcPointCheckCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcPointCheckCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcPointCheckCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_POINTCHECK_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EcPointCheckCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
