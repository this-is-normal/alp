# Reference: where the numbers come from

Commands that produce the figures each layer records.

## pcs / prover — native host

```bash
cargo run --release -p eddsa-example --features "<backend>,bb-hash,parallel" \
    --bin measure -- 253 1,4,16,64,256
```

Op counts and their in-guest cost attribution:

```bash
cargo run --release -p eddsa-example --features "<backend>,bb-hash,count-ops" \
    --bin count_verify_bf -- 253 16
```

## recursion — OpenVM guest

```bash
cargo build --release -p eddsa-example --features "<backend>,bb-hash,parallel" --bin gen_proof
./target/release/gen_proof 253 <sigs> /root/blob.bin
python3 -c 'import json;b=open("/root/blob.bin","rb").read();json.dump({"input":["0x01"+b.hex()]},open("/root/in.json","w"))'

cd <guest>
cargo openvm keygen --app-only
cargo openvm run   --mode segment --input /root/in.json          # instructions, segments; must reveal [1, …]
cargo openvm prove app --input /root/in.json --segment-max-memory $((26<<30))
cargo openvm verify app
```

`--segment-max-memory` is **not** optional in a record: OpenVM defaults to 15 GiB, which leaves
half a 32 GB GPU idle and changes the result by ~1.2×.

## aggregation

```bash
cargo openvm keygen            # app + agg_prefix
cargo openvm setup             # leaf/internal recursive pks
cargo openvm prove stark --input /root/in.json
cargo openvm verify stark
```

## Per-operation costs (for attribution, not for headlines)

`alp-openvm/guests/microbench` — second input word selects the mode; subtract mode 0.
See its README for the mode table.

## Phase attribution

`alp-openvm/guests/ablation/README.md` — the nine `ablate-*` features, the counters, and the trap:
**a cut is only valid when it removes code that exists nowhere else in the run.**
