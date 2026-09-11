#!/usr/bin/env python3
"""Basefold multilinear-PCS cost-model estimator (companion to brakedown_cost_model.py).

Two ways to set parameters:
  1) MANUAL    : edit the CONFIG dict below, run `python3 basefold_cost_model.py`.
  2) BENCHMARK : keep parameter sets in a reusable JSON file and pass it
                 `python3 basefold_cost_model.py params.json [--json]`.
                 The file overrides any CONFIG keys. `--json` emits machine-readable output.

Model (mirrors basefold/src/merkle_pcs.rs). Commit a poly of N=2^num_vars evals with a rate-1/inv_rate
foldable code (codeword_0 = inv_rate*N). Open = num_vars fold rounds: each sends a degree-2 sumcheck
univariate, folds the codeword, and Merkle-commits the folded oracle; then q queries open value+sibling
at every round oracle with Merkle paths. q = ceil((sec_bits-grind)/-log2(1-radius)) (foldable-code proximity
radius; list-decoding regime). Calibrated to measured runs (nv=12,16): proof ~0.3%, verify hashes ~0.2%,
mul/add ~1%, sqr/inv exact. Fewer calibration points than the Brakedown model — treat as first-order.
Outputs: proof size (KB/MB), #field elements, #hash codes, verifier hashes, verifier ADD/MUL/INV/SQR.
"""
import argparse
import math

# ---- CONFIG: manual defaults; a JSON file (positional arg) overrides any of these ----
CONFIG = dict(
    num_vars = 24,           # log2 of input size; N = 2**num_vars committed evaluations
    sec_bits = 100,          # target soundness in bits
    grind    = 0,            # proof-of-work bits (contributes `grind` to security -> fewer queries)
    inv_rate = 4,            # foldable-code inverse rate: codeword_0 = inv_rate*N (rate = 1/inv_rate)
    radius   = 0.215,        # proximity radius e for queries: q = ceil((sec_bits-grind)/-log2(1-e))
                             #   (list-decoding regime; ~0.215 gives q=288 at sec_bits=100 for rate 1/4)
    hash     = "babybear16",  # hash function (see HASHES)
)

# ---- Hash presets (same as the Brakedown script). Fields per hash: ----
#   perms_per_elem : Poseidon2 perms to absorb one ed25519 data element (the "hash rate")
#   digest_elems   : Merkle digest size, in hash-field elements
#   field_bytes    : bytes per hash-field element
HASHES = {
    "babybear16":     dict(perms_per_elem=9 / 8, digest_elems=8, field_bytes=4),   # SplitFieldHasher: 9 limbs / rate 8
    "ed25519_native": dict(perms_per_elem=0.5,   digest_elems=1, field_bytes=32),  # Poseidon2 width 3, rate 2
}

# ---- model constants ----
DATA_ELEM_BYTES = 32     # ed25519 field element, canonical serialized
FERMAT_SQ       = 253    # squarings in one a^(p-2) inversion (~log2 p)
PROOF_OVERHEAD  = 1.19   # serde framing of the nested per-query/per-round proof (calibrated to measured)


def num_queries(radius, sec_bits, grind=0):
    """Query count for the folding proximity test: q = ceil((sec_bits - grind) / -log2(1 - radius)).
    Returns (q, per_query_bits)."""
    per_query = -math.log2(1 - radius)               # soundness bits per query
    return max(0, math.ceil((sec_bits - grind) / per_query)), per_query


def estimate(num_vars, sec_bits, grind, inv_rate, radius, perms_per_elem, digest_elems, field_bytes):
    q, per_query = num_queries(radius, sec_bits, grind)
    d = num_vars                                     # fold rounds
    codeword0 = inv_rate << num_vars                 # base codeword length
    q = min(q, codeword0)                            # (never clamps in practice)
    achieved = q * per_query + grind
    depth0 = max(1, math.ceil(math.log2(codeword0 // 2)))     # leaves = codeword/2 (pair matrix)
    sum_depths = sum(max(1, depth0 - i) for i in range(d))    # Merkle-path total over the d oracles
    leaf_per_round = math.ceil(2 * perms_per_elem)  # a leaf is a pair (2 elements)
    grind_wit = 1 if grind > 0 else 0

    # proof: sumcheck univariates (3 each) + final value + per-query pair openings (2 per round)
    field_elems = 3 * d + 1 + q * 2 * d + grind_wit
    hash_codes = d + q * sum_depths                 # d fold roots + Merkle-path digests
    proof_bytes = round((field_elems * DATA_ELEM_BYTES + hash_codes * digest_elems * field_bytes) * PROOF_OVERHEAD)

    # verifier: replay sumcheck (no hashes) + per query: path compressions + a leaf hash per round
    v_hashes = round(q * (sum_depths + d * leaf_per_round)) + grind_wit  # +1 = PoW check
    v_inv = d                                        # one Fermat inversion per fold round (lazy twists)
    v_sqr = v_inv * FERMAT_SQ                        # each inversion ~ a^(p-2)
    v_mul = round(q * d * 3.5 + v_inv * 250)         # fold consistency + inversion muls
    v_add = q * d * 4                                # fold add/sub

    return dict(
        num_vars=num_vars, N=1 << num_vars, inv_rate=inv_rate, codeword=codeword0, fold_rounds=d,
        grind=grind, q=q, achieved_bits=round(achieved, 1), secure=(achieved >= sec_bits - 0.5),
        prover_grind_hashes=((1 << grind) if grind else 0),
        proof_bytes=proof_bytes, proof_kb=round(proof_bytes / 1000, 1), proof_mb=round(proof_bytes / 1e6, 2),
        field_elems=field_elems, hash_codes=hash_codes,
        v_hashes=v_hashes, v_mul=v_mul, v_add=v_add, v_sqr=v_sqr, v_inv=v_inv,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("params", nargs="?", help="JSON file of parameter overrides (else use CONFIG)")
    ap.add_argument("--json", action="store_true", help="machine-readable JSON output")
    a = ap.parse_args()

    p = dict(CONFIG)
    if a.params:
        import json
        with open(a.params) as fh:
            p.update(json.load(fh))

    h = HASHES[p["hash"]]
    g = lambda k, d: p.get(k, d)   # raw override if given, else hash-preset value
    res = estimate(p["num_vars"], p["sec_bits"], p["grind"], p["inv_rate"], p["radius"],
                   g("perms_per_elem", h["perms_per_elem"]), g("digest_elems", h["digest_elems"]),
                   g("field_bytes", h["field_bytes"]))

    if a.json:
        import json
        print(json.dumps(res))
    else:
        w = "" if res["secure"] else "  <-- SUB-SECURE"
        gr = f" grind={res['grind']}(prover 2^{res['grind']}={res['prover_grind_hashes']:,} hashes)" if res["grind"] else ""
        print(f"Basefold  nv={res['num_vars']} N={res['N']:,}  rate=1/{res['inv_rate']} "
              f"codeword={res['codeword']:,}  fold_rounds={res['fold_rounds']}  [{p['hash']}]{gr}")
        print(f"  queries q : {res['q']:,}   sec = {res['achieved_bits']} bits{w}")
        print(f"  proof     : {res['proof_kb']:,.1f} KB ({res['proof_mb']:,.2f} MB)   "
              f"field_elems={res['field_elems']:,}  hash_codes={res['hash_codes']:,}")
        print(f"  verifier  : hashes={res['v_hashes']:,}   MUL={res['v_mul']:,}  ADD={res['v_add']:,}"
              f"  SQR={res['v_sqr']:,}  INV={res['v_inv']:,}")


if __name__ == "__main__":
    main()
