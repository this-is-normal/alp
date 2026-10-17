//! Fibonacci AIR + trace generation over the ed25519 base field.
//!
//! Identical in shape to `multi-fib-example`'s `fib_air`, but with `F = Ed25519`.

use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_matrix::dense::RowMajorMatrix;

use crate::config::F;

/// Two columns per row: `(left, right) = (fib(i), fib(i + 1))`.
pub const NUM_COLS: usize = 2;

/// Fibonacci AIR.
///
/// - First row equals the first two public values `(a, b)`.
/// - Each transition advances the recurrence: `left' = right`, `right' = left + right`.
/// - The last row's `right` equals the claimed output `x`.
pub struct FibAir;

pub struct FibRow<T> {
    pub left: T,
    pub right: T,
}

impl<T> Borrow<FibRow<T>> for [T] {
    fn borrow(&self) -> &FibRow<T> {
        // Safety: two fields of type T in declaration order match the layout of [T; 2].
        debug_assert_eq!(self.len(), NUM_COLS);
        let ptr = self.as_ptr().cast::<FibRow<T>>();
        unsafe { &*ptr }
    }
}

impl<X> BaseAir<X> for FibAir {
    fn width(&self) -> usize {
        NUM_COLS
    }
    fn num_public_values(&self) -> usize {
        3
    }
}

impl<AB: AirBuilder> Air<AB> for FibAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let pis = builder.public_values();
        let (a, b, x) = (pis[0], pis[1], pis[2]);

        let local: &FibRow<AB::Var> = main.current_slice().borrow();
        let next: &FibRow<AB::Var> = main.next_slice().borrow();

        let mut first = builder.when_first_row();
        first.assert_eq(local.left, a);
        first.assert_eq(local.right, b);

        let mut trans = builder.when_transition();
        trans.assert_eq(local.right, next.left);
        trans.assert_eq(local.left + local.right, next.right);

        builder.when_last_row().assert_eq(local.right, x);
    }
}

/// Build a Fibonacci trace of `n` rows seeded with `(0, 1)`. `n` must be a power of two.
pub fn fib_trace(n: usize) -> RowMajorMatrix<F> {
    use p3_field::PrimeCharacteristicRing;

    let (mut left, mut right) = (F::ZERO, F::ONE);
    let mut values = Vec::with_capacity(NUM_COLS * n);
    for _ in 0..n {
        values.push(left);
        values.push(right);
        let next_left = right;
        let next_right = left + right;
        left = next_left;
        right = next_right;
    }
    RowMajorMatrix::new(values, NUM_COLS)
}

/// Public inputs `[a, b, x]` read directly from a trace (robust to field wraparound).
pub fn fib_public_values_for_trace(trace: &RowMajorMatrix<F>) -> [F; 3] {
    let n = trace.values.len() / NUM_COLS;
    let last = trace.values[(n - 1) * NUM_COLS + 1];
    [trace.values[0], trace.values[1], last]
}
