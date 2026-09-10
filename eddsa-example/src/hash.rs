//! Field, hash and transcript stack for the EdDSA STARK.
//!
//! **Hash field (feature `bb-hash`).** By default the Merkle commitment and
//! Fiat-Shamir transcript hash with a Poseidon2 permutation over **ed25519**
//! (255-bit). With `--features bb-hash` the hash instead runs over **BabyBear**
//! (31-bit): the committed data, linear code, and sumcheck stay over ed25519, and
//! only the hash field changes — each ed25519 value is decomposed into
//! base-`|BabyBear|` limbs for hashing (see [`SplitFieldHasher`]) and the
//! challenger sponges over BabyBear while presenting an ed25519 interface (see
//! [`DuplexChallengerCF`]). This makes every Poseidon2 permutation ~100x cheaper
//! inside the 32-bit OpenVM verifier; the default ed25519 hash is kept for
//! host-side comparison and op-count calibration.

use p3_ed25519::Ed25519;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};

// ---- ed25519-hash-only imports (default) ----
#[cfg(not(feature = "bb-hash"))]
use p3_challenger::DuplexChallengerPF;
#[cfg(not(feature = "bb-hash"))]
use rand::SeedableRng;
#[cfg(not(feature = "bb-hash"))]
use p3_ed25519::{
    ED25519_POSEIDON2_FULL_ROUNDS, ED25519_POSEIDON2_PARTIAL_ROUNDS_3, Poseidon2Ed25519,
};
#[cfg(not(feature = "bb-hash"))]
use rand::rngs::SmallRng;

// ---- bb-hash-only imports ----
#[cfg(feature = "bb-hash")]
use p3_baby_bear::{BabyBear, MontyBypass, Poseidon2BabyBear, default_babybear_poseidon2_16};
#[cfg(feature = "bb-hash")]
use p3_challenger::DuplexChallengerCF;
#[cfg(feature = "bb-hash")]
use p3_symmetric::SplitFieldHasher;

// ---- count-ops instrumentation (either hash field) ----
#[cfg(feature = "count-ops")]
use p3_symmetric::{CryptographicPermutation, Permutation};
/// Merkle cap height: the commitment becomes the `2^MERKLE_CAP_HEIGHT` digests
/// that many levels below the root, and every authentication path is that much
/// shorter. `MerkleTreeMmcs` clamps it per tree, so one value is safe even for
/// the shallow oracles late in a BaseFold fold chain.
///
/// **Measured: leave this at 0.** A cap looks like an easy win on paper — and it
/// is, against *unpruned* paths — but path pruning already removes the top
/// levels: with hundreds of queries the upper ~log₂(queries) levels are covered
/// by some other query and their siblings are recomputable for free. A cap then
/// buys nothing off the paths while adding `2^h` digests to **every**
/// commitment, and Quasar has many (the column tree, the block tree, the
/// coefficient tree, and one fold-oracle root per round per opening).
///
/// At cap 6 the EdDSA proof grew: 2.61 → 2.69 MB at 16 signatures, 3.33 → 3.43
/// MB at 64. Worth revisiting only if the query counts drop far enough that the
/// top levels stop being saturated.
pub(crate) const MERKLE_CAP_HEIGHT: usize = 0;

/// Base (trace) field.
pub type F = Ed25519;
/// Challenge field: the degree-1 self-extension of ed25519 (`|EF| ≈ 2^255`).
pub type EF = Ed25519;

// =====================================================================================
// Hash stack A: Poseidon2-over-ed25519 (default).
// =====================================================================================

/// Poseidon2 state width (ed25519 Poseidon2 supports width 3).
#[cfg(not(feature = "bb-hash"))]
const WIDTH: usize = 3;
/// Sponge rate (challenger + hash); capacity = 1 (~255 bits).
#[cfg(not(feature = "bb-hash"))]
const RATE: usize = 2;
/// Merkle digest / sponge output size in field elements.
#[cfg(not(feature = "bb-hash"))]
const OUT: usize = 1;

#[cfg(all(not(feature = "bb-hash"), not(feature = "count-ops")))]
type Perm = Poseidon2Ed25519<WIDTH>;
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
type Perm = CountingPerm;

/// Counts every Poseidon2-over-Ed25519 permutation (a hash invocation). Locally the hash is over
/// Ed25519, so its muls are also counted in `MUL`; counting perms lets us subtract them to isolate
/// the PCS *arithmetic* muls (with `bb-hash` the hash is BabyBear, so this subtraction is moot).
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
#[derive(Clone, Debug)]
pub struct CountingPerm(Poseidon2Ed25519<WIDTH>);
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
impl Permutation<[Ed25519; WIDTH]> for CountingPerm {
    fn permute_mut(&self, input: &mut [Ed25519; WIDTH]) {
        p3_ed25519::counters::on_perm();
        self.0.permute_mut(input);
    }
}
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
impl CryptographicPermutation<[Ed25519; WIDTH]> for CountingPerm {}

#[cfg(not(feature = "bb-hash"))]
pub(crate) type MyHash = PaddingFreeSponge<Perm, WIDTH, RATE, OUT>;
#[cfg(not(feature = "bb-hash"))]
pub(crate) type MyCompress = TruncatedPermutation<Perm, 2, OUT, WIDTH>;
#[cfg(not(feature = "bb-hash"))]
pub type MyMmcs = MerkleTreeMmcs<Ed25519, Ed25519, MyHash, MyCompress, 2, OUT>;

/// The general-`PrimeField` duplex challenger (ed25519 is not `PrimeField64`).
#[cfg(not(feature = "bb-hash"))]
pub type MyChallenger = DuplexChallengerPF<Ed25519, Perm, WIDTH, RATE>;

/// Fixed RNG seed for the (DEMO-grade) Poseidon2 round constants; seeding
/// deterministically makes the prover and verifier build the identical permutation.
#[cfg(not(feature = "bb-hash"))]
const POSEIDON2_SEED: u64 = 0xED25519_C0FFEE;

/// Build the fixed Poseidon2 permutation (ed25519 is not `PrimeField64`, so the
/// generic `new_from_rng` is used).
#[cfg(not(feature = "bb-hash"))]
fn perm() -> Perm {
    let mut rng = SmallRng::seed_from_u64(POSEIDON2_SEED);
    let inner = Poseidon2Ed25519::<WIDTH>::new_from_rng(
        ED25519_POSEIDON2_FULL_ROUNDS,
        ED25519_POSEIDON2_PARTIAL_ROUNDS_3,
        &mut rng,
    );
    #[cfg(not(feature = "count-ops"))]
    {
        inner
    }
    #[cfg(feature = "count-ops")]
    {
        CountingPerm(inner)
    }
}

/// The Merkle MMCS (leaf hash + node compression) over the ed25519 permutation.
#[cfg(not(feature = "bb-hash"))]
pub(crate) fn build_mmcs() -> MyMmcs {
    MyMmcs::new(MyHash::new(perm()), MyCompress::new(perm()), MERKLE_CAP_HEIGHT)
}

/// Calibration: ed25519 muls consumed by ONE Poseidon2 permutation (to subtract hashing from the
/// verify mul-count and isolate PCS arithmetic). Only meaningful for the ed25519 hash.
#[cfg(all(not(feature = "bb-hash"), feature = "count-ops"))]
pub fn muls_per_perm() -> u64 {
    use p3_field::PrimeCharacteristicRing;
    let p = perm();
    let mut st = [F::ZERO; WIDTH];
    p3_ed25519::counters::reset();
    p.permute_mut(&mut st);
    p3_ed25519::counters::snapshot().mul
}

/// A bare seeded challenger — identical on the prove and verify sides (no domain separator).
#[cfg(not(feature = "bb-hash"))]
pub fn challenger() -> MyChallenger {
    MyChallenger::new(perm())
}

// =====================================================================================
// Hash stack B: Poseidon2-over-BabyBear (feature `bb-hash`).
// =====================================================================================

/// BabyBear Poseidon2 state width / sponge rate / digest size (31-bit hash field).
#[cfg(feature = "bb-hash")]
const BB_WIDTH: usize = 16;
#[cfg(feature = "bb-hash")]
const BB_RATE: usize = 8;
/// Merkle digest = 8 BabyBear ≈ 248-bit (~124-bit collision resistance).
#[cfg(feature = "bb-hash")]
const BB_DIGEST: usize = 8;

/// The permutation is `MontyBypass`-wrapped: `P'(v) = R⁻¹·P(R·v)`.
///
/// Cryptographically it is `P` conjugated by a field bijection — same strength.
/// Operationally it is what lets an OpenVM guest feed the Poseidon2 precompile
/// the raw Montgomery words with **no per-lane conversion**, which measured 549
/// instructions per permutation against 16 for the bare chip call. The host
/// pays two elementwise constant multiplications, which stay vectorised.
#[cfg(all(feature = "bb-hash", not(feature = "count-ops")))]
type PermBB = MontyBypass<Poseidon2BabyBear<BB_WIDTH>>;
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
type PermBB = CountingPerm;

/// Wraps the BabyBear Poseidon2 permutation to count every invocation (hash permutation)
/// during profiling. Value-identical to the inner permutation.
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
#[derive(Clone, Debug)]
pub struct CountingPerm(Poseidon2BabyBear<BB_WIDTH>);
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
impl Permutation<[BabyBear; BB_WIDTH]> for CountingPerm {
    fn permute_mut(&self, input: &mut [BabyBear; BB_WIDTH]) {
        p3_ed25519::counters::on_hash();
        self.0.permute_mut(input);
    }
}
#[cfg(all(feature = "bb-hash", feature = "count-ops"))]
impl CryptographicPermutation<[BabyBear; BB_WIDTH]> for CountingPerm {}

/// Inner BabyBear→BabyBear sponge (hashes the limb stream).
#[cfg(feature = "bb-hash")]
type InnerHash = PaddingFreeSponge<PermBB, BB_WIDTH, BB_RATE, BB_DIGEST>;
/// Leaf hasher: ed25519 elements → BabyBear digest via base-|BabyBear| limbs.
#[cfg(feature = "bb-hash")]
pub(crate) type MyHash = SplitFieldHasher<Ed25519, BabyBear, InnerHash>;
#[cfg(feature = "bb-hash")]
pub(crate) type MyCompress = TruncatedPermutation<PermBB, 2, BB_DIGEST, BB_WIDTH>;
#[cfg(feature = "bb-hash")]
pub type MyMmcs = MerkleTreeMmcs<Ed25519, BabyBear, MyHash, MyCompress, 2, BB_DIGEST>;

/// Transcript challenger: sponges over BabyBear, observes/samples ed25519.
#[cfg(feature = "bb-hash")]
pub type MyChallenger = DuplexChallengerCF<Ed25519, BabyBear, PermBB, BB_WIDTH, BB_RATE>;

/// Build the fixed BabyBear Poseidon2 permutation.
#[cfg(feature = "bb-hash")]
fn perm_bb() -> PermBB {
    // OpenVM-matching Horizen RC16 constants (identical to the poseidon2 precompile).
    #[cfg(not(feature = "count-ops"))]
    {
        MontyBypass::new(default_babybear_poseidon2_16())
    }
    #[cfg(feature = "count-ops")]
    {
        CountingPerm(MontyBypass::new(default_babybear_poseidon2_16()))
    }
}

/// The Merkle MMCS: SplitFieldHasher leaves + BabyBear-permutation node compression.
#[cfg(feature = "bb-hash")]
pub(crate) fn build_mmcs() -> MyMmcs {
    MyMmcs::new(
        MyHash::new(InnerHash::new(perm_bb())),
        MyCompress::new(perm_bb()),
        MERKLE_CAP_HEIGHT,
    )
}

/// A bare seeded challenger — identical on the prove and verify sides (no domain separator).
#[cfg(feature = "bb-hash")]
pub fn challenger() -> MyChallenger {
    MyChallenger::new(perm_bb())
}

/// The stacked-witness layout the multi-stark opening path uses (prefix binding
/// order). Shared by every commitment backend.
pub type L = p3_sumcheck::layout::PrefixProver<F, EF>;
