use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Trace columns for one twisted-Edwards point addition `(x1,y1) + (x2,y2) = (x3,y3)`.
///
/// `p`, `q`, `m` are helper witness columns that keep every constraint at degree <= 2:
/// - `p = x1*x2`
/// - `q = y1*y2`
/// - `m = d*p*q`
///
/// with `a = -1` fixed by the curve, the complete addition law becomes:
/// - `x3*(1+m) = x1*y2 + y1*x2`
/// - `y3*(1-m) = q + p`
///
/// This same layout and formula covers doubling (`(x1,y1) == (x2,y2)`) with no case split.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcAddCols<T> {
    pub x1: T,
    pub y1: T,
    pub x2: T,
    pub y2: T,
    pub x3: T,
    pub y3: T,
    pub p: T,
    pub q: T,
    pub m: T,
}

pub const NUM_EC_ADD_COLS: usize = size_of::<EcAddCols<u8>>();

impl<T> Borrow<EcAddCols<T>> for [T] {
    fn borrow(&self) -> &EcAddCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_ADD_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EcAddCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcAddCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcAddCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_ADD_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EcAddCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
