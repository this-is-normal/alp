# Full pipeline at 256 signatures: 256 ed25519 sigs -> a 534,636 B aggregated proof in ~34 s of proving

<!-- layer: pipeline (prover + recursion + aggregation) -->
<!-- target: Quasar end-to-end, native prove -> prove app -> prove stark -> prove evm -->

## Target

| | |
|---|---|
| layer | `pipeline` — spans `prover`, `recursion` and `aggregation` |
| subject | one batch of 256 ed25519 signatures carried through every proving stage |
| **parameters** | `num_bits=253`, 256 signatures, trace `65536x45`, `c=2`, `BF_INV_RATE=4`, `POW_BITS=16`, `MERKLE_CAP_HEIGHT=0`, `MAX_LOG_COLS=16` (256 fits the baked table) |
| build features | `quasar,bb-hash,parallel` (host) · `rv32i,rv32m,io,modular,poseidon2` chips (guest) |
| commands | stage 1 `cargo run --release -p eddsa-example --features "quasar,bb-hash,parallel" --bin measure -- 253 256` · stages 2-4 `cargo openvm prove {app,stark,evm} --exe <vmexe> --input /root/qz_in_256.json --segment-max-memory $((26<<30))` |
| input artefacts | `/root/eddsa_quasar_256.bin` (2,214,994 B) hex-wrapped to `qz_in_256.json`; guest `.vmexe` 93 MB, `exe commit 0x005f1052…4a3a` |

## Provenance

| field | value |
|---|---|
| date | 2026-09-11 |
| tree | `alp2` @ `70cb077d` (`quazar2`) |
| **uncommitted** | **true (19 files)** — the direct-preprocessed work; not reproducible from the commit alone |
| host | vultr4 |
| CPU | AMD EPYC 7B12 — 256 visible, quota 30.72 cores |
| memory | 503 GB |
| GPU | NVIDIA GeForce RTX 5090, 32607 MiB, 580.159.03 |
| rustc | 1.98.1 |
| cargo-openvm | cargo-openvm v2.0.2 (59a69b8) [cuda] |

## Results

| stage | what it proves | wall-clock (vultr4) | output bytes |
|---|---|---|---|
| 1 — native Quasar prove | 256 ed25519 signatures | setup 254.9 ms · **prove 3,741.4 ms** · verify 124.1 ms | 2,196,130 (8,578.6 B/sig) |
| 2 — `prove app` | an OpenVM guest verified stage 1's proof | **28 s** | 7,028,421 |
| 3 — `prove stark` | recursive aggregation of stage 2 | **30 s** (cumulative, see Notes) | **534,636** |
| 4 — `prove evm` | halo2/KZG wrapper | **blocked** — see below | — |

All figures **measured**.

**Recursion inflates before it compresses.** Verifying a 2.20 MB proof inside the VM costs a
7.03 MB app proof — 3.2x *larger* than what it verifies. Only aggregation pays that back, at
534,636 B. Stage 2 in isolation is a proof-size regression and is worth doing only as a step
toward stage 3 or 4.

**Stage 2 reproduced exactly**: 28 s and 7,028,421 B, byte-identical to the same measurement in
`2026-09-11-recursion-quasar-direct-preprocessed`, from a separate invocation.

### Stage 4 is blocked on infrastructure, not on the scheme

`cargo openvm setup --evm --download` fetches the KZG SRS from S3. It completed the non-EVM
artifacts (`internal_recursive.{pk,vk}`, `root.pk` — which is why stage 3 ran) and then stalled:

```
Downloading kzg_bn254_19.srs
Error: failed to download kzg_bn254_19.srs body from S3
  0: streaming error
  1: minimum throughput was specified at 1 B/s, but throughput of 0 B/s was observed
real 13m37.941s
```

64 MB of ~2.15 GB fetched (files are `2^k · 64 B`, k=10..24), 130 GB disk free, so this is the
box's network, not capacity. A retry loop is in progress; `.srs` files already present are skipped,
so attempts resume rather than restart. **This is a vast.ai result, not a Quasar result.**

## Comparability

- **comparable to:** `2026-09-11-recursion-quasar-direct-preprocessed` — same box, commit,
  parameters and `.vmexe`; stage 2 here is that record's 256-signature row, independently repeated.
- **not comparable to:** `2026-09-07-pcs-native-comparison` for stage 1 — different box (vultr3)
  **and** before the preprocessed-commitment removal. Both differences favour this run; the 6.96 s
  -> 3.74 s gap must not be attributed to either cause alone.

## Tried and rejected

| change | result | why it lost |
|---|---|---|
| running the SRS download concurrently with stages 1-3 | kept | download is network-bound at ~0-200 KB/s and ~0% CPU; no measurable contention. Disclosed because "measured while something else ran" is a field this format exists to capture |
| waiting for the download before running any stage | reverted | would have idled the box for hours before stage 1 |

## Notes

**Stages 2-4 are cumulative, not additive.** `prove stark` and `prove evm` do not consume the
previous stage's proof file — each takes the exe and input and redoes the chain internally. So
stage 3's 30 s *contains* stage 2's 28 s. Never sum this column.

**Stage 3's ~2 s marginal cost over stage 2 now REPRODUCES at a second batch size.** At 256
signatures `prove app` is 28 s and `prove stark` 30 s; at 16 signatures they are 21 s and 23 s.
The same +2 s at two sizes, with the rest scaling as the app proof does, is consistent with
`prove stark` redoing the app work plus a small constant GPU aggregation. It is still a
surprisingly cheap aggregation and worth a closer look, but it is repeatable, not an artefact.

**The root proof is 534,636 B, not the 630,732 B documented in `alp-openvm/README.md` and
`guests/README.md`. RESOLVED: the "constant" claim is right; the value is stale.** Tested by
proving at a second batch size — 16 signatures gives **534,636 B**, byte-identical to 256's, in
23 s. So the root proof is genuinely independent of batch size, as documented; 630,732 B simply
belongs to the 2026-09-09 epoch (different guest, and `internal_recursive.pk`/`root.pk` have since
been rewritten by the `--evm` setup attempts). Which of those moved the constant is still not
isolated, and cannot be from here: the old keys were overwritten in place and the 2026-09-09
`stark.proof` was overwritten by stage 3. **The docs should quote 534,636 B for this toolchain.**

**Stage 1's script-level `secs=40` is almost entirely `cargo` compilation** and is not reported
above; the `measure` binary's own setup/prove/verify split is the prover characterisation.

**Two driver bugs cost time here, both the same mistake:** matching a prefix true of success *and*
failure. `grep -q SETUPDONE` matched `SETUPDONE_1` (a failure), so stage 4 launched against absent
artefacts; then a monitor's stop condition matched that stage-4 failure line and stopped early.
Gate on the exit code, not on the presence of a marker.
