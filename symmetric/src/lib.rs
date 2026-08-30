//! A framework for symmetric cryptography primitives.

#![no_std]

extern crate alloc;

mod compression;
mod hash;
mod hasher;
mod permutation;
mod serializing_hasher;
mod split_hash;
mod sponge;

pub use compression::*;
pub use hash::*;
pub use hasher::*;
pub use permutation::*;
pub use serializing_hasher::*;
pub use split_hash::*;
pub use sponge::*;
