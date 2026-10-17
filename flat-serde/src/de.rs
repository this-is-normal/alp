//! The reader. Mirrors [`crate::ser`] exactly; change one and you change both.
//!
//! The point of the whole crate lives in [`Deserializer::deserialize_bytes`]:
//! it hands the visitor a slice **borrowed from the input**, so a 32-byte field
//! element costs a bounds check and a pointer, not a copy into a fresh `Vec`.

use serde::Deserialize;
use serde::de::{self, DeserializeSeed, IntoDeserializer, SeqAccess, Visitor};

use crate::error::{Error, Result};

/// Deserialize a value that borrows from `input`.
///
/// # Errors
///
/// If the input is truncated, malformed for `T`, or has bytes left over.
pub fn from_bytes<'de, T: Deserialize<'de>>(input: &'de [u8]) -> Result<T> {
    let mut de = Deserializer { input };
    let value = T::deserialize(&mut de)?;
    if de.input.is_empty() {
        Ok(value)
    } else {
        Err(Error::TrailingBytes(de.input.len()))
    }
}

/// Reads the flat encoding from a byte slice, borrowing where it can.
pub struct Deserializer<'de> {
    input: &'de [u8],
}

impl<'de> Deserializer<'de> {
    /// Wrap a slice. Use [`from_bytes`] unless you need partial reads.
    #[must_use]
    pub const fn new(input: &'de [u8]) -> Self {
        Self { input }
    }

    /// Bytes not yet consumed.
    #[must_use]
    pub const fn remaining(&self) -> &'de [u8] {
        self.input
    }

    #[inline]
    fn take(&mut self, n: usize) -> Result<&'de [u8]> {
        if self.input.len() < n {
            return Err(Error::UnexpectedEnd);
        }
        let (head, tail) = self.input.split_at(n);
        self.input = tail;
        Ok(head)
    }

    #[inline]
    fn len(&mut self) -> Result<usize> {
        let bytes = self.take(4)?;
        Ok(u32::from_le_bytes(bytes.try_into().expect("four bytes")) as usize)
    }
}

/// Fixed little-endian reads, matching `ser`'s `primitive!`.
macro_rules! primitive {
    ($($method:ident, $ty:ty, $visit:ident;)*) => {$(
        #[inline]
        fn $method<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
            const N: usize = core::mem::size_of::<$ty>();
            let bytes = self.take(N)?;
            let value = <$ty>::from_le_bytes(bytes.try_into().expect("fixed width"));
            visitor.$visit(value)
        }
    )*};
}

impl<'de> de::Deserializer<'de> for &mut Deserializer<'de> {
    type Error = Error;

    primitive! {
        deserialize_i8, i8, visit_i8; deserialize_i16, i16, visit_i16;
        deserialize_i32, i32, visit_i32; deserialize_i64, i64, visit_i64;
        deserialize_i128, i128, visit_i128;
        deserialize_u8, u8, visit_u8; deserialize_u16, u16, visit_u16;
        deserialize_u32, u32, visit_u32; deserialize_u64, u64, visit_u64;
        deserialize_u128, u128, visit_u128;
        deserialize_f32, f32, visit_f32; deserialize_f64, f64, visit_f64;
    }

    /// Not self-describing by design — the schema is the Rust type.
    fn deserialize_any<V: Visitor<'de>>(self, _visitor: V) -> Result<V::Value> {
        Err(Error::NotSelfDescribing)
    }

    fn deserialize_bool<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        match self.take(1)?[0] {
            0 => visitor.visit_bool(false),
            1 => visitor.visit_bool(true),
            _ => Err(Error::InvalidTag),
        }
    }

    fn deserialize_char<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        let bytes = self.take(4)?;
        let raw = u32::from_le_bytes(bytes.try_into().expect("four bytes"));
        visitor.visit_char(char::from_u32(raw).ok_or(Error::OutOfRange)?)
    }

    fn deserialize_str<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        let n = self.len()?;
        let bytes = self.take(n)?;
        visitor.visit_borrowed_str(core::str::from_utf8(bytes).map_err(|_| Error::InvalidUtf8)?)
    }

    fn deserialize_string<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        self.deserialize_str(visitor)
    }

    /// **The hot path.** Borrowed, so a field element is parsed in place.
    #[inline]
    fn deserialize_bytes<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        let n = self.len()?;
        let bytes = self.take(n)?;
        visitor.visit_borrowed_bytes(bytes)
    }

    fn deserialize_byte_buf<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        self.deserialize_bytes(visitor)
    }

    fn deserialize_option<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        match self.take(1)?[0] {
            0 => visitor.visit_none(),
            1 => visitor.visit_some(self),
            _ => Err(Error::InvalidTag),
        }
    }

    fn deserialize_unit<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        visitor.visit_unit()
    }

    fn deserialize_unit_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value> {
        visitor.visit_unit()
    }

    fn deserialize_newtype_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        visitor: V,
    ) -> Result<V::Value> {
        visitor.visit_newtype_struct(self)
    }

    fn deserialize_seq<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        let n = self.len()?;
        visitor.visit_seq(Elements { de: self, remaining: n })
    }

    /// Tuples and fixed arrays carry no length — both sides know it.
    fn deserialize_tuple<V: Visitor<'de>>(self, len: usize, visitor: V) -> Result<V::Value> {
        visitor.visit_seq(Elements { de: self, remaining: len })
    }

    fn deserialize_tuple_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        len: usize,
        visitor: V,
    ) -> Result<V::Value> {
        self.deserialize_tuple(len, visitor)
    }

    fn deserialize_map<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        let n = self.len()?;
        visitor.visit_map(Elements { de: self, remaining: n })
    }

    fn deserialize_struct<V: Visitor<'de>>(
        self,
        _name: &'static str,
        fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value> {
        self.deserialize_tuple(fields.len(), visitor)
    }

    fn deserialize_enum<V: Visitor<'de>>(
        self,
        _name: &'static str,
        _variants: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value> {
        visitor.visit_enum(self)
    }

    fn deserialize_identifier<V: Visitor<'de>>(self, visitor: V) -> Result<V::Value> {
        self.deserialize_u32(visitor)
    }

    fn deserialize_ignored_any<V: Visitor<'de>>(self, _visitor: V) -> Result<V::Value> {
        Err(Error::NotSelfDescribing)
    }

    fn is_human_readable(&self) -> bool {
        false
    }
}

/// A run of values: sequence elements, tuple fields, or map entries.
struct Elements<'a, 'de> {
    de: &'a mut Deserializer<'de>,
    remaining: usize,
}

impl<'de> SeqAccess<'de> for Elements<'_, 'de> {
    type Error = Error;

    #[inline]
    fn next_element_seed<T: DeserializeSeed<'de>>(&mut self, seed: T) -> Result<Option<T::Value>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        self.remaining -= 1;
        seed.deserialize(&mut *self.de).map(Some)
    }

    fn size_hint(&self) -> Option<usize> {
        Some(self.remaining)
    }
}

impl<'de> de::MapAccess<'de> for Elements<'_, 'de> {
    type Error = Error;

    fn next_key_seed<K: DeserializeSeed<'de>>(&mut self, seed: K) -> Result<Option<K::Value>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        self.remaining -= 1;
        seed.deserialize(&mut *self.de).map(Some)
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(&mut self, seed: V) -> Result<V::Value> {
        seed.deserialize(&mut *self.de)
    }

    fn size_hint(&self) -> Option<usize> {
        Some(self.remaining)
    }
}

impl<'de> de::EnumAccess<'de> for &mut Deserializer<'de> {
    type Error = Error;
    type Variant = Self;

    fn variant_seed<V: DeserializeSeed<'de>>(self, seed: V) -> Result<(V::Value, Self)> {
        let bytes = self.take(4)?;
        let index = u32::from_le_bytes(bytes.try_into().expect("four bytes"));
        let value = seed.deserialize(index.into_deserializer())?;
        Ok((value, self))
    }
}

impl<'de> de::VariantAccess<'de> for &mut Deserializer<'de> {
    type Error = Error;

    fn unit_variant(self) -> Result<()> {
        Ok(())
    }

    fn newtype_variant_seed<T: DeserializeSeed<'de>>(self, seed: T) -> Result<T::Value> {
        seed.deserialize(self)
    }

    fn tuple_variant<V: Visitor<'de>>(self, len: usize, visitor: V) -> Result<V::Value> {
        de::Deserializer::deserialize_tuple(self, len, visitor)
    }

    fn struct_variant<V: Visitor<'de>>(
        self,
        fields: &'static [&'static str],
        visitor: V,
    ) -> Result<V::Value> {
        de::Deserializer::deserialize_tuple(self, fields.len(), visitor)
    }
}
