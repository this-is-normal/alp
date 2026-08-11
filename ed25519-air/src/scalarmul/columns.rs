use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

use crate::add::NUM_EC_ADD_COLS;

/// Trace columns for one double-and-add step of scalar multiplication, `row_i` corresponding
/// to bit `i` of the scalar (LSB first).
///
/// Invariants maintained across the trace (see [`crate::scalarmul::air::EcScalarMulAir`]):
/// - `(x_dbl, y_dbl)` holds `2^i * P` at row `i`.
/// - `(x_acc, y_acc)` holds `sum(bit_k * 2^k * P for k in 0..=i-1)` at row `i` (the running
///   result, *before* this row's bit is folded in).
/// - `weight` holds `2^i`.
/// - `scalar_acc` holds `sum(bit_k * 2^k for k in 0..=i)` (cumulative, *including* this row's bit).
///
/// The `add_*` and `dbl_*` blocks are laid out to exactly match [`crate::add::EcAddCols`] so
/// they can be checked by embedding [`crate::add::EcAddAir`] as a sub-AIR twice per row: once
/// for `acc + dbl` (the candidate next accumulator, folded in only if `bit = 1`) and once for
/// `dbl + dbl` (doubling, always folded in).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcScalarMulCols<T> {
    pub bit: T,

    // -- add sub-air block: x1,y1,x2,y2 = (x_acc,y_acc,x_dbl,y_dbl); x3,y3 = candidate acc+dbl --
    pub x_acc: T,
    pub y_acc: T,
    pub x_dbl: T,
    pub y_dbl: T,
    pub sum_x: T,
    pub sum_y: T,
    pub add_p: T,
    pub add_q: T,
    pub add_m: T,

    // -- double sub-air block: x1,y1,x2,y2 = copies of (x_dbl,y_dbl); x3,y3 = 2*dbl --
    pub dbl_x1: T,
    pub dbl_y1: T,
    pub dbl_x2: T,
    pub dbl_y2: T,
    pub dbl_x3: T,
    pub dbl_y3: T,
    pub dbl_p: T,
    pub dbl_q: T,
    pub dbl_m: T,

    /// `2^i`, doubles every row.
    pub weight: T,
    /// Running weighted sum of bits seen so far (including this row's bit).
    pub scalar_acc: T,

    /// Accumulator value after folding in this row's bit (`mux(bit, sum, acc)`); equals
    /// `next.x_acc`/`next.y_acc` for every row but the last, where it's the only witnessed
    /// copy of the final `scalar * base` result (there's no row after the last to read it
    /// from otherwise).
    pub result_x: T,
    pub result_y: T,
}

pub const NUM_EC_SCALARMUL_COLS: usize = size_of::<EcScalarMulCols<u8>>();

// The two sub-air blocks above must exactly match `EcAddCols`'s layout/width.
const _: () = assert!(9 == NUM_EC_ADD_COLS);

impl<T> Borrow<EcScalarMulCols<T>> for [T] {
    fn borrow(&self) -> &EcScalarMulCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_SCALARMUL_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EcScalarMulCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcScalarMulCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcScalarMulCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_SCALARMUL_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EcScalarMulCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

/// Column index (within the flat row) where the "add" sub-air block starts.
pub const ADD_BLOCK_START: usize = 1;
/// Column index (within the flat row) where the "double" sub-air block starts.
pub const DBL_BLOCK_START: usize = ADD_BLOCK_START + NUM_EC_ADD_COLS;
