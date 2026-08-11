use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

use crate::add::{EcAddCols, NUM_EC_ADD_COLS};
use crate::fixed_scalarmul::{EcFixedScalarMulCols, NUM_EC_FIXED_SCALARMUL_COLS};
use crate::scalarmul::{EcScalarMulCols, NUM_EC_SCALARMUL_COLS};
use crate::subgroupcheck::{EcSubgroupCheckCols, NUM_EC_SUBGROUPCHECK_COLS};

/// Width of the *region*: the leading columns of every main-trace row, reinterpreted according
/// to the row type selected by the preprocessed trace.
///
/// | row type | region contents |
/// |---|---|
/// | ladder | [`EdDsaLadderCols`] (both scalar-mul ladders side by side) |
/// | subgroup check | [`EcSubgroupCheckCols`] |
/// | add | [`EcAddCols`], remaining columns unused |
/// | padding | unused |
///
/// Sized by the widest occupant. The two ladders happen to pack to exactly the same width as
/// one subgroup check, so nothing is wasted on the rows that matter.
pub const EDDSA_REGION_WIDTH: usize = 35;

const _: () = assert!(NUM_EDDSA_LADDER_COLS <= EDDSA_REGION_WIDTH);
const _: () = assert!(NUM_EC_SUBGROUPCHECK_COLS <= EDDSA_REGION_WIDTH);
const _: () = assert!(NUM_EC_ADD_COLS <= EDDSA_REGION_WIDTH);

/// Region layout on a ladder row: the two double-and-add ladders run in lockstep, one scalar bit
/// per row, sharing the row but nothing else.
///
/// - `var` computes `[k]A` against a witnessed base (tied to the subgroup-checked `A` on row 0).
/// - `fixed` computes `[s]B` against the generator, reading `2^i * B` from the preprocessed
///   trace instead of witnessing a doubling block.
///
/// Both sub-layouts are reused verbatim from the standalone AIRs so the constraint logic can be
/// transcribed field-for-field; see [`crate::eddsa::air::EdDsaAir`].
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EdDsaLadderCols<T> {
    pub var: EcScalarMulCols<T>,
    pub fixed: EcFixedScalarMulCols<T>,
}

pub const NUM_EDDSA_LADDER_COLS: usize = size_of::<EdDsaLadderCols<u8>>();

// The nested blocks must exactly match the standalone AIRs' layouts/widths.
const _: () = assert!(NUM_EDDSA_LADDER_COLS == NUM_EC_SCALARMUL_COLS + NUM_EC_FIXED_SCALARMUL_COLS);

/// Column index (within the flat row) where the variable-base ladder block starts.
pub const VAR_LADDER_START: usize = 0;
/// Column index (within the flat row) where the fixed-base ladder block starts.
pub const FIXED_LADDER_START: usize = VAR_LADDER_START + NUM_EC_SCALARMUL_COLS;

/// Values that must agree across non-adjacent rows. An AIR only ever sees `local` and `next`, so
/// anything wired between, say, the subgroup-check row for `A` and ladder row 0 has to ride along
/// on every row in between. These are held constant by a degree-1 `next.v == local.v` transition
/// constraint each, and are the sole mechanism tying the row types together.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EdDsaSharedCols<T> {
    /// Public key `A`, subgroup-checked on its own row and used as the variable ladder's base.
    pub a_x: T,
    pub a_y: T,
    /// Signature component `R`, subgroup-checked on its own row and used as the add row's `P1`.
    pub r_x: T,
    pub r_y: T,
    /// Signature component `s`, reconstructed from the fixed ladder's bits.
    pub s: T,
    /// Challenge scalar `k`, reconstructed from the variable ladder's bits.
    pub k: T,
    /// `[s]B`, produced by the fixed ladder's last row and consumed by the final equality.
    pub sb_x: T,
    pub sb_y: T,
    /// `[k]A`, produced by the variable ladder's last row and consumed as the add row's `P2`.
    pub ka_x: T,
    pub ka_y: T,
}

pub const NUM_EDDSA_SHARED_COLS: usize = size_of::<EdDsaSharedCols<u8>>();

/// Column index (within the flat row) where the shared columns start.
pub const SHARED_START: usize = EDDSA_REGION_WIDTH;

/// Total main-trace width: the reinterpreted region followed by the shared columns.
pub const NUM_EDDSA_COLS: usize = EDDSA_REGION_WIDTH + NUM_EDDSA_SHARED_COLS;

/// Public values, in trace order: constrained on the first row against the shared columns of the
/// same name. `[s]B`, `[k]A` and the witness points are *not* public.
pub const NUM_EDDSA_PUBLIC_VALUES: usize = 6;

/// Preprocessed columns: fixed at setup, identical for every proof against a given
/// `(d, base_b, num_bits)`.
///
/// The selectors are what make the single table polymorphic — every constraint in
/// [`crate::eddsa::air::EdDsaAir`] is multiplied by exactly one of them, so a row type's
/// constraints vanish on every other row type (and on padding, where all selectors are zero).
/// Because they live here rather than in the main trace, they need no constraints of their own:
/// a prover cannot choose them.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct EdDsaPreprocessedCols<T> {
    /// `2^i * B` on ladder row `i`, feeding the fixed-base ladder's add law. Junk (the identity)
    /// on non-ladder rows, where no constraint reads it.
    pub dbl_x: T,
    pub dbl_y: T,

    /// Rows `0 .. num_bits`: both ladders step here.
    pub is_ladder: T,
    /// Rows `0 .. num_bits - 1`: the ladder rows that have a valid *next* ladder row. A dedicated
    /// column rather than `is_ladder * next.is_ladder` so transition gating stays a single
    /// degree-1 factor.
    pub is_ladder_step: T,
    /// Row `num_bits - 1`: the ladder's last row, where the results are read out.
    pub is_ladder_last: T,
    /// Row `num_bits`: subgroup check on `A`.
    pub is_sgc_a: T,
    /// Row `num_bits + 1`: subgroup check on `R`.
    pub is_sgc_r: T,
    /// Row `num_bits + 2`: `R + [k]A`, and the final equality against `[s]B`.
    pub is_add: T,
}

pub const NUM_EDDSA_PREPROCESSED_COLS: usize = size_of::<EdDsaPreprocessedCols<u8>>();

/// The two-ladder view of a row's region.
pub fn ladder_cols<T>(row: &[T]) -> &EdDsaLadderCols<T> {
    row[..NUM_EDDSA_LADDER_COLS].borrow()
}

pub fn ladder_cols_mut<T>(row: &mut [T]) -> &mut EdDsaLadderCols<T> {
    row[..NUM_EDDSA_LADDER_COLS].borrow_mut()
}

/// The subgroup-check view of a row's region.
pub fn sgc_cols<T>(row: &[T]) -> &EcSubgroupCheckCols<T> {
    row[..NUM_EC_SUBGROUPCHECK_COLS].borrow()
}

pub fn sgc_cols_mut<T>(row: &mut [T]) -> &mut EcSubgroupCheckCols<T> {
    row[..NUM_EC_SUBGROUPCHECK_COLS].borrow_mut()
}

/// The point-addition view of a row's region.
pub fn add_cols<T>(row: &[T]) -> &EcAddCols<T> {
    row[..NUM_EC_ADD_COLS].borrow()
}

pub fn add_cols_mut<T>(row: &mut [T]) -> &mut EcAddCols<T> {
    row[..NUM_EC_ADD_COLS].borrow_mut()
}

/// The shared columns, present on *every* row regardless of type.
pub fn shared_cols<T>(row: &[T]) -> &EdDsaSharedCols<T> {
    row[SHARED_START..NUM_EDDSA_COLS].borrow()
}

pub fn shared_cols_mut<T>(row: &mut [T]) -> &mut EdDsaSharedCols<T> {
    row[SHARED_START..NUM_EDDSA_COLS].borrow_mut()
}

impl<T> Borrow<EdDsaLadderCols<T>> for [T] {
    fn borrow(&self) -> &EdDsaLadderCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_LADDER_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EdDsaLadderCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EdDsaLadderCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EdDsaLadderCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_LADDER_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EdDsaLadderCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

impl<T> Borrow<EdDsaSharedCols<T>> for [T] {
    fn borrow(&self) -> &EdDsaSharedCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_SHARED_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EdDsaSharedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EdDsaSharedCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EdDsaSharedCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_SHARED_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EdDsaSharedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

impl<T> Borrow<EdDsaPreprocessedCols<T>> for [T] {
    fn borrow(&self) -> &EdDsaPreprocessedCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_PREPROCESSED_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<EdDsaPreprocessedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<EdDsaPreprocessedCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut EdDsaPreprocessedCols<T> {
        debug_assert_eq!(self.len(), NUM_EDDSA_PREPROCESSED_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<EdDsaPreprocessedCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
