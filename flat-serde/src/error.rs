//! One error type for both directions.

use alloc::string::{String, ToString};
use core::fmt::{self, Display};

/// What can go wrong reading or writing the flat format.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// The input ended mid-value.
    UnexpectedEnd,
    /// A length or discriminant did not fit its fixed width.
    OutOfRange,
    /// A `bool` byte was neither 0 nor 1, or an `Option` tag was invalid.
    InvalidTag,
    /// Invalid UTF-8 in a `str`.
    InvalidUtf8,
    /// The format is not self-describing, so `deserialize_any` cannot work.
    NotSelfDescribing,
    /// Bytes remained after the value was fully read.
    TrailingBytes(usize),
    /// A message from serde itself.
    Message(String),
}

impl Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEnd => f.write_str("input ended mid-value"),
            Self::OutOfRange => f.write_str("length or discriminant out of range"),
            Self::InvalidTag => f.write_str("invalid bool or option tag"),
            Self::InvalidUtf8 => f.write_str("invalid utf-8"),
            Self::NotSelfDescribing => {
                f.write_str("the flat format is not self-describing; deserialize_any is unsupported")
            }
            Self::TrailingBytes(n) => write!(f, "{n} trailing bytes"),
            Self::Message(m) => f.write_str(m),
        }
    }
}

impl core::error::Error for Error {}

impl serde::ser::Error for Error {
    fn custom<T: Display>(msg: T) -> Self {
        Self::Message(msg.to_string())
    }
}

impl serde::de::Error for Error {
    fn custom<T: Display>(msg: T) -> Self {
        Self::Message(msg.to_string())
    }
}

/// Shorthand for this crate's results.
pub type Result<T> = core::result::Result<T, Error>;
