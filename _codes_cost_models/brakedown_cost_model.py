#!/usr/bin/env python3
"""Brakedown-ML PCS cost-model estimator.

Two ways to set parameters:
  1) MANUAL    : edit the CONFIG dict below, run `python3 brakedown_cost_model.py`.
  2) BENCHMARK : keep parameter sets in a reusable JSON file and pass it
                 `python3 brakedown_cost_model.py params.json [--json]`.
                 The file overrides any CONFIG keys (num_vars, sec_bits, conjectured, grind, split,
                 preset, hash, and/or raw rho/delta/tn_over_n/perms_per_elem/digest_elems/field_bytes).
                 `--json` emits machine-readable output so a benchmark harness can diff these
                 theoretical numbers against its own practical results after each run.

Code/hash pick a preset (PRESETS / HASHES); the JSON can also override any numeric param directly.
Model: see BENCHMARKS.md. Outputs: proof size (KB/MB), #field elements, #hash codes, verifier
hashes, verifier ADD/MUL/INV/SQR.
"""
import argparse
import math

# ---- CONFIG: manual defaults; each key is overridable by the CLI flag of the same name ----
CONFIG = dict(
    num_vars    = 24,             # log2 of input size; N = 2**num_vars committed evaluations
    sec_bits    = 100,            # target soundness in bits
    conjectured = False,          # False = provable radius e=delta/3 ; True = conjectured e=delta/2 (DP Conj-1)
    grind       = 0,              # proof-of-work bits; contributes `grind`, so t only needs (sec_bits - grind)
    split       = "min_size",     # dimension strategy (see STRATEGIES / choose_log_rows):
                                  #   square_matrix | paper | min_size | min_arith | min_hashes | min_verify
    preset      = "small-proof",  # code operating point (see PRESETS)
    hash        = "babybear16",   # hash function (see HASHES)
    hash_cost   = 1.0,            # min_verify weight: guest cost per verifier hash (perm)
    op_cost     = 1.0,            # min_verify weight: guest cost per verifier field op (MUL/ADD)
)

# ---- Code presets (Haböck Table 1 / lcpc SdigCode). Fields per preset: ----
#   alpha     : recursion split  (code construction; documentation only)
#   rho       : code rate = message/codeword ; inverse rate r = 1/rho          [used by model]
#   delta     : relative minimum distance = beta/r                             [used by model]
#   c_a, c_b  : sparse-block row weights c_n / d_n  (documentation only)
#   tn_over_n : linear-encode field-ops per element (verifier re-encode cost)  [used by model]
PRESETS = {
    "small-proof": dict(alpha=0.238, rho=0.581, delta=0.070, c_a=9, c_b=12, tn_over_n=18.3),
    "balanced":    dict(alpha=0.200, rho=0.610, delta=0.050, c_a=7, c_b=12, tn_over_n=13.7),
    "fast-prover": dict(alpha=0.120, rho=0.704, delta=0.020, c_a=6, c_b=16, tn_over_n=9.9),
}

# ---- Hash presets. Fields per hash: ----
#   perms_per_elem : Poseidon2 perms to absorb one ed25519 data element (the "hash rate")
#   digest_elems   : Merkle digest size, in hash-field elements
#   field_bytes    : bytes per hash-field element
HASHES = {
    "babybear16":     dict(perms_per_elem=9 / 8, digest_elems=8, field_bytes=4),   # SplitFieldHasher: 9 limbs / rate 8
    "ed25519_native": dict(perms_per_elem=0.5,   digest_elems=1, field_bytes=32),  # Poseidon2 width 3, rate 2
}

# ---- model constants ----
DATA_ELEM_BYTES = 32     # ed25519 field element, canonical serialized
COMBOS          = 1      # combined rows in the proof: 1 after the Diamond–Posen test/eval consolidation
PROOF_OVERHEAD  = 1.14   # serde framing beyond the element count; calibrated to the consolidated+optimal
                         #   impl (~±5% across sizes; set 1.0 for the raw lower bound)


STRATEGIES = ("square_matrix", "paper", "min_size", "min_arith", "min_hashes", "min_verify")


def split_metrics(num_vars, lr, t, rho, tn_over_n, perms_per_elem, digest_elems, field_bytes):
    """Cost metrics of one row/column split (log_rows = lr) at query count t."""
    rows = 1 << lr
    cols = 1 << (num_vars - lr)
    codeword = math.ceil(cols / rho)                    # encoded columns = Merkle leaves
    depth = max(1, math.ceil(math.log2(codeword)))
    leaf = math.ceil(rows * perms_per_elem)             # per-query leaf hash (a rows-tall column)
    field_elems = t * rows + COMBOS * cols
    hash_codes = t * depth
    proof_bytes = round((field_elems * DATA_ELEM_BYTES + hash_codes * digest_elems * field_bytes) * PROOF_OVERHEAD)
    v_hashes = round(t * (leaf + depth) + COMBOS * cols * perms_per_elem)
    v_arith = 2 * (COMBOS * t * rows + round((tn_over_n / 2 + 1) * cols))  # MUL + ADD
    return dict(rows=rows, cols=cols, codeword=codeword, depth=depth,
                field_elems=field_elems, hash_codes=hash_codes,
                proof_bytes=proof_bytes, v_hashes=v_hashes, v_arith=v_arith)


def choose_log_rows(num_vars, t, rho, tn_over_n, perms_per_elem, digest_elems, field_bytes, strategy,
                    hash_cost=1.0, op_cost=1.0):
    """Pick log_rows for the tensor split. All strategies keep codeword_len >= t (enough columns
    to run t queries → the target security is reachable):
      square_matrix : log_rows = num_vars//2 (equal-ish rows/cols).
      paper         : closed-form balanced-proof rule rows ≈ sqrt(N/t)  (Ligero/Brakedown), i.e.
                      log_rows = round((num_vars - log2 t)/2), widened until codeword >= t.
      min_size      : brute-force minimise the actual proof BYTES.
      min_arith     : brute-force minimise verifier field ops (MUL + ADD).
      min_hashes    : brute-force minimise verifier hash (Poseidon2 perm) count.
      min_verify    : brute-force minimise weighted guest verify cost
                      hash_cost*v_hashes + op_cost*v_arith (the recursion-relevant total).
    """
    def enough(lr):
        return math.ceil((1 << (num_vars - lr)) / rho) >= t

    # square_matrix is the explicit naive baseline: it may be insecure (estimate() reports
    # secure=False for it) — that is the point of keeping it for comparison. Every other
    # strategy TARGETS the security level, so if no shape can run t queries it is an error.
    unreachable = (f"{strategy}: no split of a 2^{num_vars} poly can run {t} queries at rate "
                   f"rho={rho:.3g} — even the widest codeword has < {t} columns; "
                   f"increase num_vars or pick a higher-delta preset")
    if strategy == "square_matrix":
        return num_vars // 2
    if strategy == "paper":
        lr = max(0, min(num_vars, round((num_vars - math.log2(max(t, 1))) / 2)))
        while lr > 0 and not enough(lr):
            lr -= 1                                       # widen (more columns) until enough
        if not enough(lr):
            raise ValueError(unreachable)
        return lr
    costfn = {
        "min_size":   lambda m: m["proof_bytes"],
        "min_arith":  lambda m: m["v_arith"],
        "min_hashes": lambda m: m["v_hashes"],
        "min_verify": lambda m: hash_cost * m["v_hashes"] + op_cost * m["v_arith"],
    }[strategy]
    best_lr, best = None, None
    for lr in range(num_vars + 1):
        if not enough(lr):
            continue
        c = costfn(split_metrics(num_vars, lr, t, rho, tn_over_n, perms_per_elem, digest_elems, field_bytes))
        if best is None or c < best:
            best, best_lr = c, lr
    if best_lr is None:
        raise ValueError(unreachable)                     # assert, never silently return insecure
    return best_lr


def num_queries(delta, sec_bits, conjectured, grind=0):
    """Column-query count for the tensor proximity test:
        t = ceil( (sec_bits - grind) / -log2(1 - e) ),  e = delta/3 (provable) or delta/2 (conjectured).
    Grinding supplies `grind` bits of security for free, so fewer queries are needed.
    Returns (t, e, per_query_bits)."""
    e = delta / 2 if conjectured else delta / 3       # tested proximity radius
    per_query = -math.log2(1 - e)                     # soundness bits per opened column
    return max(0, math.ceil((sec_bits - grind) / per_query)), e, per_query


def estimate(num_vars, sec_bits, conjectured, grind, split, rho, delta, tn_over_n,
             perms_per_elem, digest_elems, field_bytes, hash_cost=1.0, op_cost=1.0):
    t_req, e, per_query = num_queries(delta, sec_bits, conjectured, grind)
    lr = choose_log_rows(num_vars, t_req, rho, tn_over_n, perms_per_elem, digest_elems, field_bytes, split,
                         hash_cost, op_cost)
    codeword_cols = math.ceil((1 << (num_vars - lr)) / rho)
    t = min(t_req, codeword_cols)                      # queries, clamped to available columns
    m = split_metrics(num_vars, lr, t, rho, tn_over_n, perms_per_elem, digest_elems, field_bytes)
    achieved = t * per_query + grind                    # actual bits achieved
    grind_wit = 1 if grind > 0 else 0                  # PoW nonce carried in the proof
    proof_bytes = m["proof_bytes"] + round(grind_wit * DATA_ELEM_BYTES * PROOF_OVERHEAD)
    v_muladd = m["v_arith"] // 2                        # per side (MUL == ADD in this model)
    return dict(
        num_vars=num_vars, N=1 << num_vars, split=split, rows=m["rows"], cols=m["cols"], codeword=m["codeword"],
        regime=("d/2" if conjectured else "d/3"), grind=grind, t=t, achieved_bits=round(achieved, 1),
        secure=(achieved >= sec_bits - 0.5), prover_grind_hashes=((1 << grind) if grind else 0),
        proof_bytes=proof_bytes, proof_kb=round(proof_bytes / 1000, 1), proof_mb=round(proof_bytes / 1e6, 2),
        field_elems=m["field_elems"] + grind_wit, hash_codes=m["hash_codes"],
        v_hashes=m["v_hashes"] + grind_wit, v_mul=v_muladd, v_add=v_muladd, v_sqr=0, v_inv=0,  # no sqr/inv
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

    c, h = PRESETS[p["preset"]], HASHES[p["hash"]]
    g = lambda k, d: p.get(k, d)   # raw override if given, else preset value
    res = estimate(p["num_vars"], p["sec_bits"], p["conjectured"], p["grind"], p["split"],
                   g("rho", c["rho"]), g("delta", c["delta"]), g("tn_over_n", c["tn_over_n"]),
                   g("perms_per_elem", h["perms_per_elem"]), g("digest_elems", h["digest_elems"]),
                   g("field_bytes", h["field_bytes"]), p["hash_cost"], p["op_cost"])

    if a.json:
        import json
        print(json.dumps(res))
    else:
        w = "" if res["secure"] else "  <-- SUB-SECURE (t clamped to codeword)"
        gr = f" grind={res['grind']}(prover 2^{res['grind']}={res['prover_grind_hashes']:,} hashes)" if res["grind"] else ""
        print(f"Brakedown-ML  nv={res['num_vars']} N={res['N']:,}  {res['split']} split "
              f"{res['rows']:,}x{res['cols']:,}  codeword={res['codeword']:,}  "
              f"[{p['preset']}/{p['hash']}/{res['regime']}]{gr}")
        print(f"  queries t : {res['t']:,}   sec = {res['achieved_bits']} bits{w}")
        print(f"  proof     : {res['proof_kb']:,.1f} KB ({res['proof_mb']:,.2f} MB)   "
              f"field_elems={res['field_elems']:,}  hash_codes={res['hash_codes']:,}")
        print(f"  verifier  : hashes={res['v_hashes']:,}   MUL={res['v_mul']:,}  ADD={res['v_add']:,}"
              f"  SQR={res['v_sqr']}  INV={res['v_inv']}")


if __name__ == "__main__":
    main()
