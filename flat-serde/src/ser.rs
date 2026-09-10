//! The writer. Mirrors [`crate::de`] exactly; change one and you change both.

use alloc::vec::Vec;

use serde::{Serialize, ser};

use crate::error::{Error, Result};

/// Serialize to a fresh buffer.
///
/// # Errors
///
/// Only from `T`'s own `Serialize`, or a length exceeding `u32`.
pub fn to_allocvec<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>> {
    let mut ser = Serializer { out: Vec::new() };
    value.serialize(&mut ser)?;
    Ok(ser.out)
}

/// Writes the flat encoding into a growing buffer.
pub struct Serializer {
    out: Vec<u8>,
}

impl Serializer {
    fn push_len(&mut self, len: usize) -> Result<()> {
        let len = u32::try_from(len).map_err(|_| Error::OutOfRange)?;
        self.out.extend_from_slice(&len.to_le_bytes());
        Ok(())
    }
}

/// Fixed little-endian, no varints.
macro_rules! primitive {
    ($($method:ident, $ty:ty;)*) => {$(
        fn $method(self, v: $ty) -> Result<()> {
            self.out.extend_from_slice(&v.to_le_bytes());
            Ok(())
        }
    )*};
}

impl ser::Serializer for &mut Serializer {
    type Ok = ();
    type Error = Error;
    type SerializeSeq = Self;
    type SerializeTuple = Self;
    type SerializeTupleStruct = Self;
    type SerializeTupleVariant = Self;
    type SerializeMap = Self;
    type SerializeStruct = Self;
    type SerializeStructVariant = Self;

    primitive! {
        serialize_i8, i8; serialize_i16, i16; serialize_i32, i32; serialize_i64, i64;
        serialize_u8, u8; serialize_u16, u16; serialize_u32, u32; serialize_u64, u64;
        serialize_f32, f32; serialize_f64, f64;
    }

    fn serialize_bool(self, v: bool) -> Result<()> {
        self.out.push(u8::from(v));
        Ok(())
    }

    fn serialize_i128(self, v: i128) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }

    fn serialize_u128(self, v: u128) -> Result<()> {
        self.out.extend_from_slice(&v.to_le_bytes());
        Ok(())
    }

    /// `usize` travels as `u64` so a 32-bit guest and a 64-bit host agree.
    fn serialize_char(self, v: char) -> Result<()> {
        self.serialize_u32(v as u32)
    }

    fn serialize_str(self, v: &str) -> Result<()> {
        self.serialize_bytes(v.as_bytes())
    }

    fn serialize_bytes(self, v: &[u8]) -> Result<()> {
        self.push_len(v.len())?;
        self.out.extend_from_slice(v);
        Ok(())
    }

    fn serialize_none(self) -> Result<()> {
        self.out.push(0);
        Ok(())
    }

    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<()> {
        self.out.push(1);
        value.serialize(self)
    }

    fn serialize_unit(self) -> Result<()> {
        Ok(())
    }

    fn serialize_unit_struct(self, _name: &'static str) -> Result<()> {
        Ok(())
    }

    fn serialize_unit_variant(
        self,
        _name: &'static str,
        index: u32,
        _variant: &'static str,
    ) -> Result<()> {
        self.serialize_u32(index)
    }

    fn serialize_newtype_struct<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        value: &T,
    ) -> Result<()> {
        value.serialize(self)
    }

    fn serialize_newtype_variant<T: Serialize + ?Sized>(
        self,
        _name: &'static str,
        index: u32,
        _variant: &'static str,
        value: &T,
    ) -> Result<()> {
        self.serialize_u32(index)?;
        value.serialize(self)
    }

    fn serialize_seq(self, len: Option<usize>) -> Result<Self::SerializeSeq> {
        self.push_len(len.ok_or(Error::OutOfRange)?)?;
        Ok(self)
    }

    /// Tuples and fixed arrays have a length known to both sides — no prefix.
    fn serialize_tuple(self, _len: usize) -> Result<Self::SerializeTuple> {
        Ok(self)
    }

    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleStruct> {
        Ok(self)
    }

    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleVariant> {
        self.serialize_u32(index)?;
        Ok(self)
    }

    fn serialize_map(self, len: Option<usize>) -> Result<Self::SerializeMap> {
        self.push_len(len.ok_or(Error::OutOfRange)?)?;
        Ok(self)
    }

    fn serialize_struct(self, _name: &'static str, _len: usize) -> Result<Self::SerializeStruct> {
        Ok(self)
    }

    fn serialize_struct_variant(
        self,
        _name: &'static str,
        index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeStructVariant> {
        self.serialize_u32(index)?;
        Ok(self)
    }

    fn is_human_readable(&self) -> bool {
        false
    }
}

/// Every compound is just its elements back to back.
macro_rules! compound {
    ($($trait:ident, $method:ident $(, $key:ident)?;)*) => {$(
        impl ser::$trait for &mut Serializer {
            type Ok = ();
            type Error = Error;

            $(
                fn $key<T: Serialize + ?Sized>(&mut self, key: &T) -> Result<()> {
                    key.serialize(&mut **self)
                }
            )?

            fn $method<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<()> {
                value.serialize(&mut **self)
            }

            fn end(self) -> Result<()> {
                Ok(())
            }
        }
    )*};
}

compound! {
    SerializeSeq, serialize_element;
    SerializeTuple, serialize_element;
    SerializeTupleStruct, serialize_field;
    SerializeTupleVariant, serialize_field;
    SerializeMap, serialize_value, serialize_key;
}

/// Struct fields carry no names on the wire.
macro_rules! struct_compound {
    ($($trait:ident;)*) => {$(
        impl ser::$trait for &mut Serializer {
            type Ok = ();
            type Error = Error;

            fn serialize_field<T: Serialize + ?Sized>(
                &mut self,
                _key: &'static str,
                value: &T,
            ) -> Result<()> {
                value.serialize(&mut **self)
            }

            fn end(self) -> Result<()> {
                Ok(())
            }
        }
    )*};
}

struct_compound! { SerializeStruct; SerializeStructVariant; }
