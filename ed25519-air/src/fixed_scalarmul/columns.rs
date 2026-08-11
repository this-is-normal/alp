use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Trace columns for one double-and-add step of scalar multiplication against a fixed base
/// point, `row_i` corresponding to bit `i` of the scalar (LSB first).
///
/// Invariants maintained across the trace (see
/// [`crate::fixed_scalarmul::air::EcFixedScalarMulAir`]):
/// - `(x_acc, y_acc)` holds `sum(bit_k * 2^k * base for k in 0..=i-1)` at row `i` (the running
///   result, *before* this row's bit is folded in).
/// - `weight` holds `2^i`.
/// - `scalar_acc` holds `sum(bit_k * 2^k for k in 0..=i)` (cumulative, *including* this row's bit).
///
/// Unlike [`crate::scalarmul::columns::EcScalarMulCols`], there is no doubling-point witness
/// block: `2^i * base` is read straight out of [`EcFixedScalarMulPreprocessedCols`] instead of
/// being re-derived and proved on-circuit every row.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcFixedScalarMulCols<T> {
    pub bit: T,

    // -- add-law block: (x_acc,y_acc) + (x_dbl,y_dbl, from preprocessed) = (sum_x,sum_y) --
    pub x_acc: T,
    pub y_acc: T,
    pub sum_x: T,
    pub sum_y: T,
    pub add_p: T,
    pub add_q: T,
    pub add_m: T,

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

pub const NUM_EC_FIXED_SCALARMUL_COLS: usize = size_of::<EcFixedScalarMulCols<u8>>();

impl<T> Borrow<EcFixedScalarMulCols<T>> for [T] {
    fn borrow(&self) -> &EcFixedScalarMulCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_FIXED_SCALARMUL_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EcFixedScalarMulCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcFixedScalarMulCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcFixedScalarMulCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_FIXED_SCALARMUL_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EcFixedScalarMulCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

/// Preprocessed columns: `2^i * base` at row `i`, generated once per `(base, num_rows)` and
/// shared across every proof against that base (not committed per-proof).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EcFixedScalarMulPreprocessedCols<T> {
    pub x_dbl: T,
    pub y_dbl: T,
}

pub const NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS: usize =
    size_of::<EcFixedScalarMulPreprocessedCols<u8>>();

impl<T> Borrow<EcFixedScalarMulPreprocessedCols<T>> for [T] {
    fn borrow(&self) -> &EcFixedScalarMulPreprocessedCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS);
        let (prefix, shorts, suffix) =
            unsafe { self.align_to::<EcFixedScalarMulPreprocessedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EcFixedScalarMulPreprocessedCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EcFixedScalarMulPreprocessedCols<T> {
        debug_assert_eq!(self.len(), NUM_EC_FIXED_SCALARMUL_PREPROCESSED_COLS);
        let (prefix, shorts, suffix) =
            unsafe { self.align_to_mut::<EcFixedScalarMulPreprocessedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
