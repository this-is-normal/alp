# Changes and optimisations

What this fork adds on top of Plonky3, clustered by area. Each entry is the change, why, and the
measured effect. Reverted experiments are kept at the bottom — they cost real time to discover.

Companion log for the zkVM side: `code/alp-openvm/OPTIMIZATIONS.md` (chips and guests).
Individual benchmark runs with full provenance: [`bench_results/`](bench_results/).

Numbers are **measured** unless marked *(modelled)*. Where a figure is "of the guest", it means the
OpenVM recursion verifier at 16 signatures.

---

## PCS / Quasar — a new scheme

**New crates `qa-code/` and `quasar/`** — quasi-abelian code commitment plus code switching into
BaseFold (eprint 2026/1839). WHT-based encoder (additions only, no twiddles), five sumchecks, two
batched BaseFold openings.
*Fastest native prover in the tree: 6.96 s at 256 signatures against Brakedown's 7.44 s, with a
3.2× smaller proof and 4.0× faster verification.*

**`security/src/qa.rs`** — the [LLX⁺26b] distance bound and the error terms. Reproduces the paper's
Table 4 exactly, all ten rows.
*Establishes that QA codes need big fields: δ→0 at 31 bits, 0.46 at 255 — ed25519 is the best case.*

**`QuasarPcs::new_verifier` / `QuasarMlPcs::new_verifier`** — construct a verifier from the
coefficient *commitment* instead of the QA code. `new` ends with a full encode-and-Merkle-commit
that the verifier only reads a digest from.
*Removes an O(N₂) setup per proof.*

**`backend::quasar::resolved_params`** — the host resolves `QuasarParams::choose` and ships the
shape in the proof blob. The chooser bisects on `f64` `log`/`pow`/`expm1`.
*−95 M instructions: a circuit with no FPU was running softfloat to rediscover four integers.*

## PCS / BaseFold

**Batched openings** (`commit_many`, `open_batch`, `verify_batch`) — one opening covering several
polynomials, which is what code switching needs.

**Montgomery batch inversion in `RandomFoldableCode::new`** — the twist vector was inverted one
element at a time; over ed25519 an inversion is 252 squarings plus 77 multiplications.
*Setup 22.2 s → 0.81 s at 16 signatures.*

**Addressable sampling + `new_verifier` / `diag_at`** — each code element derives from its own
ChaCha seed, so a verifier can reach position `(level, j)` without replaying the stream.
*Makes an O(q·log N) verifier possible instead of O(N) — and lets the parameters be baked into a
guest image entirely.*

**`new_verifier_with_tables`** — take twists **and their inverses** from a caller-supplied table.
*Without the inverses the lazy path pays a 14,879-instruction field inversion per queried position,
which is what made it lose to eager construction twice.*

**16-bit grinding** (`new_with_grinding`) and a **parallel deterministic search**
(`grind_parallel`) returning the smallest witness, so proofs stay bit-identical across thread
counts.
*−14 % proof size; the parallel search took Quasar's prove from 9.77 s to 6.96 s, because six
serial grinds were ~55 % of the prover.*

## PCS / Brakedown and Ligero

**Merkle path pruning** through `Mmcs::{open_multi_batch, verify_multi_batch}` — queries share
ancestors, so the shared siblings are recomputable and need not be sent.
*−61 % proof and −55 % verify for Quasar at 256 signatures; helps in proportion to how
path-dominated a proof is.*

**Ligero backend** — Reed–Solomon in the tensor PCS: high distance, few queries, but superlinear
encoding.

## PCS / shared Merkle

**`verify_batch` fast path** for one matrix, arity 2, cap 0, power-of-two height — skips the sort,
the allocated arity schedule and the injection logic.
*20,295 → 14,958 instructions per path in a microbenchmark. No effect on the recursion guests,
which use the pruned path; kept for other callers.*

## Field / ed25519

**`from_canonical_bytes` / `to_canonical_bytes`** replacing the serde implementation. The old one
allocated a `Vec<u8>`, built two `BigUint`s, **rebuilt the modulus as a `BigUint` on every call**,
and round-tripped through bytes again.
*1,830 → ~60 instructions per element; −92 M of the guest, with no wire-format change.*

**`canonical_u64_words`** (new `PrimeField` method, overridden here) — the canonical value as four
`u64`s with no allocation.
*Feeds the limb split below; on a zkVM guest the canonical form is already the internal one.*

**`StandardUniform` without `BigUint`** — sampling built a `BigUint` from 32 bytes and unpacked it
again.
*5,127 → 4,008 instructions per element; every foldable-code twist costs one.*

**`#[cfg_attr(target_os = "zkvm", repr(C, align(32)))]`** — OpenVM's modular precompile operates on
a 32-byte-aligned type, so every operation *copied* `[u64; 4]` in and out. Matching the alignment
makes the conversion a pointer cast.
*Multiply 178 → 66, square 186 → **34** (the bare chip op). −5.4 M of the guest — far less than the
op-counts predicted.*

**`from_canonical_words`** — const construction from the canonical form, free on a guest.
*Lets a baked parameter table be read for ~34 instructions instead of 6,809 to re-derive.*

## Field / BabyBear and monty-31

**`MontyParameters::MONTY_R2` and `to_monty` via Montgomery multiply** — `to_monty` did
`u64 % u64`, which a 32-bit target resolves to a software `__umoddi3` call.
*Part of taking the limb split from 1,583 to 222 instructions.*

**`MontyField31::{to_monty_raw, from_monty_raw}`** — representation-level access to the Montgomery
word.

**`MontyBypass`** (`baby-bear/src/monty_bypass.rs`) — defines `P'(v) = R⁻¹·P(R·v)`. Since `R·v` is
the field element whose *canonical* value is `v`'s Montgomery word, a guest computes `P'` by handing
the raw words to the Poseidon2 precompile with **no conversion**. `P` conjugated by a bijection, so
identical strength.
*A wrapped permutation cost 549 instructions against 16 for the bare chip call. −12 % of the guest.
Changes the hash, so host and guest must move together.*

## Hashing / cross-field

**`split_pf_to_packed_limbs_into`** — the ed25519→BabyBear limb split writes into a caller stack
buffer, and gets its words from `canonical_u64_words` instead of shifting a `BigUint` per limb.
*1,583 → 222 instructions per hashed element — every Merkle leaf and transcript observation pays
one.*

## Serialisation

**New crate `flat-serde/` (`p3-flat-serde`)** — fixed-width integers, `u32` length prefixes, no
self-description, and `deserialize_bytes` handing the visitor **borrowed** bytes so a field element
is parsed in place. Keeps all the derived structure, so no proof type changed.
*postcard decodes a `Vec<Ed25519>` at 710 instructions per element against 152 flat. −17.4 M of the
guest; the blob grows ~5 %.*

**`VerifyingKey: Serialize + Deserialize`** — so a recursion guest receives the preprocessed
commitment instead of regenerating and re-committing the preprocessed trace in-circuit.

## multi-stark

**`EdDsaConfig::pcs_pair`** — read out public setup material (Quasar's coefficient commitments) a
guest would otherwise recompute.

**`config_for_eddsa_with_security`** — explicit soundness target, since Brakedown's query count is
now distance-derived and 100 bits needs more columns than a small test trace has.

## Parallelism

**rayon throughout** `basefold` and `quasar`: encode recursion (`join` above level 8), fold,
`eq_table`, `bind_msb`, the fused per-round `h(0)/h(1)/h(2)` pass, pair-matrix, tensor table,
`combine_rows`, `claim_weights`.
*Quasar 8.1× on ~23 cores. Gotcha: maybe-rayon's two backends disagree on `reduce`, so use a
`Sums<F, N>` newtype with `Add`+`Sum` and call `.sum()`.*

## Instrumentation

**Op counters** — `p3_ed25519::counters` (mul/square/add/inverse/perm), `p3_field::split_counter`,
`p3_merkle_tree::pruned_counter`, `p3_quasar::opening_probe`. All behind `count-ops`.

**Ablation features** — nine `ablate-*` flags across `multi-stark`, `basefold` and `quasar` that
make `verify` return early at a phase boundary, so a guest's cost splits by phase.
*Cost measurement only — they accept unconditionally. One trap: a cut is only valid when it removes
code that exists nowhere else; cutting between two calls to the same function let the compiler
specialise the survivor and produced a phantom 3.4× asymmetry.*

---

## Tried and rejected

Kept so they are not re-attempted.

| change | result | why it lost |
|---|---|---|
| Merkle caps (height 6) | proof 2.61 → 2.69 MB | pruning already frees the top levels; a cap adds 2^h digests to **every** commitment |
| Quasar inner rate 1/4 → 1/16 | 1.26 B → 2.13 B guest instructions | 4× longer inner codeword ⇒ 4× the eager twist sampling; also slows the native prover, which is Quasar's advantage |
| lazy foldable code for Quasar (without baked inverses) | 557 M → 822 M | a field inversion per queried position beat the sampling it saved; Quasar's inner code is small enough that `q·log N₂` exceeds `N₂` |
| allocation-free streaming `LimbIter` in `SplitFieldHasher` | path 20,264 → 20,789; bulk hashing +12 % | the guest's bump allocator is nearly free; the per-item state machine is not |
| baking twists without their inverses | 177 M → 368 M | sampling became free, the per-position inversion remained |
