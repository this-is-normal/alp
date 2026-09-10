//! A flat, fixed-width serde format for data a zkVM guest has to parse.
//!
//! # Why
//!
//! A recursion guest deserialises its own input before it can verify anything,
//! and that input is almost entirely field elements. Measured in an OpenVM
//! guest, postcard decodes a `Vec<Ed25519>` at **710 instructions per element**
//! against **152** for a hand-rolled flat read — and the gap is not the integer
//! encoding, it is per-element machinery: varint decoding, bounds plumbing and
//! error paths on every primitive.
//!
//! This format removes all of it:
//!
//! * integers are fixed little-endian (`usize` travels as `u64`),
//! * lengths are a `u32` prefix, never a varint,
//! * byte strings are handed to the visitor **borrowed**, so a 32-byte field
//!   element is parsed straight out of the input with no copy,
//! * nothing is self-describing — the schema is the Rust type on both sides.
//!
//! It is deliberately not a general-purpose format. It is not portable across
//! type changes, has no versioning, and `deserialize_any` is an error. Both
//! sides must be built from the same types.
//!
//! # Trade
//!
//! Encodings are slightly larger than postcard's (fixed-width lengths, no
//! varint compaction). For proof blobs that is noise: the payload is 32-byte
//! field elements, which no format compresses.

#![no_std]

extern crate alloc;

mod de;
mod error;
mod ser;

pub use de::{Deserializer, from_bytes};
pub use error::{Error, Result};
pub use ser::{Serializer, to_allocvec};
