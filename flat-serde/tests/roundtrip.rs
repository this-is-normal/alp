//! The format must round-trip everything a proof blob contains, and must
//! reject truncation rather than reading past the end.

use p3_ed25519::Ed25519;
use p3_field::PrimeCharacteristicRing;
use p3_flat_serde::{Error, from_bytes, to_allocvec};
use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Nested {
    flag: bool,
    count: u64,
    name: String,
    rows: Vec<Vec<Ed25519>>,
    fixed: [Ed25519; 3],
    maybe: Option<Ed25519>,
    nothing: Option<u32>,
    tagged: Tag,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
enum Tag {
    Unit,
    Newtype(u64),
    Tuple(u32, Ed25519),
    Struct { a: u8, b: Vec<u16> },
}

fn element(i: u64) -> Ed25519 {
    Ed25519::from_u64(i) + Ed25519::ONE
}

fn sample(tag: Tag) -> Nested {
    Nested {
        flag: true,
        count: u64::MAX,
        name: "quasar".into(),
        rows: vec![vec![element(1), element(2)], vec![], vec![element(3)]],
        fixed: [element(4), element(5), element(6)],
        maybe: Some(element(7)),
        nothing: None,
        tagged: tag,
    }
}

#[test]
fn round_trips_every_shape_a_proof_uses() {
    for tag in [
        Tag::Unit,
        Tag::Newtype(9),
        Tag::Tuple(3, element(11)),
        Tag::Struct { a: 7, b: vec![1, 2, 3] },
    ] {
        let value = sample(tag);
        let bytes = to_allocvec(&value).expect("serialize");
        let back: Nested = from_bytes(&bytes).expect("deserialize");
        assert_eq!(value, back);
    }
}

/// Field elements are the payload, so `Vec<Ed25519>` gets its own check at a
/// size where a per-element bug would show.
#[test]
fn round_trips_a_large_field_element_vector() {
    let v: Vec<Ed25519> = (0..4096).map(element).collect();
    let bytes = to_allocvec(&v).expect("serialize");
    // 4 bytes of outer length, then a 4-byte length and 32 bytes per element.
    assert_eq!(bytes.len(), 4 + 4096 * (4 + 32));
    let back: Vec<Ed25519> = from_bytes(&bytes).expect("deserialize");
    assert_eq!(v, back);
}

#[test]
fn truncation_is_an_error_not_a_panic() {
    let v: Vec<Ed25519> = (0..8).map(element).collect();
    let bytes = to_allocvec(&v).expect("serialize");
    for cut in [0, 1, 4, 20, bytes.len() - 1] {
        let err = from_bytes::<Vec<Ed25519>>(&bytes[..cut]).unwrap_err();
        assert!(
            matches!(err, Error::UnexpectedEnd | Error::Message(_)),
            "cut at {cut} gave {err:?}"
        );
    }
}

#[test]
fn trailing_bytes_are_rejected() {
    let mut bytes = to_allocvec(&element(1)).expect("serialize");
    bytes.push(0);
    assert!(matches!(
        from_bytes::<Ed25519>(&bytes).unwrap_err(),
        Error::TrailingBytes(1)
    ));
}

/// A non-canonical field element must still be refused — the encoding is
/// looser than postcard's, the validation is not.
#[test]
fn non_canonical_field_elements_are_rejected() {
    let mut bytes = to_allocvec(&element(1)).expect("serialize");
    for b in &mut bytes[4..] {
        *b = 0xff;
    }
    assert!(from_bytes::<Ed25519>(&bytes).is_err());
}
