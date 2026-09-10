#!/usr/bin/env python3
"""Quasar PCS cost-model estimator (companion to basefold_cost_model.py / brakedown_cost_model.py).

Two ways to set parameters:
  1) MANUAL    : edit the CONFIG dict below, run `python3 quasar_cost_model.py`.
  2) BENCHMARK : `python3 quasar_cost_model.py params.json [--json]`, where the file is either a
                 dict of CONFIG overrides or a LIST of such dicts. `--json` emits machine-readable
                 output so a harness can diff it against measured numbers.

Validation: `python3 quasar_cost_model.py --check bench_data/quasar_measured.json` diffs every
prediction against the profiler's exact counters. Regenerate the measurements with
`cargo run -p pcs-profile --release --bin quasar_estimate -- 14:6:2,16:6:2,...` (each config is
num_vars:log_rows:index). As of 2026-09-06, over 16 configurations spanning num_vars 12-20,
log_rows 4-8 and index 2/4:

    t, hash_leaf, hash_compress, hash_sponge, inv, sq, mul   EXACT
    add                                                      <= 0.30%
    proof_bytes                                              <= 0.99%

Two independent halves:

  * SECURITY  — the QA-code distance bound of [LLX+26b] (Thm 4.4 / Cor 3.21) and the query count and
    union bound Quasar builds on it. Mirrors `security/src/qa.rs` op for op, and reproduces the
    paper's Table 4 EXACTLY (delta_QA, t, lambda_dist/samp/fld/total, all ten rows).

  * COST      — a structural model of the implementation in `quasar/src/pcs.rs`, for the ed25519 data
    + BabyBear-Poseidon2 hash configuration that `pcs-profile` measures. The combinatorial
    quantities (hash permutations by role, field inversions, squarings, proof structure) are derived
    from the code, not fitted: they match the measured counters EXACTLY. Verifier mul/add carry two
    fitted coefficients for per-inversion and per-round arithmetic (documented at MUL_PER_INV etc.);
    proof bytes carry one for postcard's varint digest encoding.

Prover cost is quoted as operation counts only (commitment mults/adds/hashes per Table 1 of the
paper); wall-clock is not modelled -- measure it.
"""
import argparse
import json
import math
import sys

# ---- CONFIG: manual defaults; a JSON file overrides any of these ----
CONFIG = dict(
    num_vars     = 20,     # log2 of input size; N = 2**num_vars committed evaluations
    log_rows     = 6,      # n1: matrix rows N1 = 2**log_rows (paper's own benchmarks use 6)
    index        = 2,      # QA index c = inverse rate; codeword row = c * N2
    field_bits   = 255,    # 255 = ed25519 (this fork); 127 = the paper's Mersenne127
    sec_bits     = 103,    # target bits for the column-sampling error (Table 4 uses 103)
    dist_failure = 103,    # target bits for the QA distance bound's own failure probability
    bf_inv_rate  = 4,      # BaseFold inverse rate for the code-switching openings
    bf_queries   = 288,    # BaseFold query repetitions (list-decoding regime, ~100 bits at rate 1/4)
    hash         = "babybear16",  # hash preset (see HASHES)
)

# ---- Hash presets. Fields per hash: ----
#   limbs        : small-field limbs one ed25519 element splits into (SplitFieldHasher / challenger)
#   rate         : sponge rate, in small-field elements
#   digest_elems : Merkle digest size, in small-field elements
#   digest_bytes : postcard-serialized bytes per digest (see DIGEST_BYTES_NOTE)
HASHES = {
    # SplitFieldHasher<Ed25519, BabyBear> + PaddingFreeSponge<_, 16, 8, 8>, the pcs-profile stack.
    "babybear16": dict(limbs=9, rate=8, digest_elems=8, digest_bytes=39.0),
    # Poseidon2 over ed25519 directly (width 3, rate 2, 1-element digest), the standalone config.
    "ed25519_native": dict(limbs=1, rate=2, digest_elems=1, digest_bytes=33.0),
}

# ---- model constants ----
# All of these are MEASURED, not fitted: `cargo run -p pcs-profile --release --bin
# ed25519_primitives` prints them. Only ADD_CONST is a residual fit (see below).
#
# An ed25519 element serializes as `serialize_bytes(&[u8; 32])` -> postcard writes a 1-byte varint
# length then 32 bytes.
FIELD_BYTES = 33
# DIGEST_BYTES: a BabyBear digest limb is a Montgomery u32 under postcard's varint, so its size is
# data-dependent (5 bytes for ~87% of uniform values, 4 for most of the rest). 39 B over 8 limbs is
# the measured mean; it is the one place the byte model is statistical rather than exact.
VEC_FRAMING = 2.0
# Fermat inversion a^(p-2) over ed25519.
SQ_PER_INV = 252
INV_MUL = 77
# One `RoundPolyInterpolator`: construction, then per round. `eval` runs TWO batch inversions.
INTERP_NEW = {2: (92, 6), 3: (113, 12)}       # degree -> (mul, add)
INTERP_EVAL = {2: (170, 10), 3: (187, 13)}    # degree -> (mul, add), 2 inversions each
# `basefold::interp_deg2`: a fresh `TWO.try_inverse()` every round (77 mul + 252 sq) plus 8 muls.
BF_INTERP_MUL = INV_MUL + 8
BF_INTERP_ADD = 6
# Sampling one ed25519 challenge recombines `limbs` small-field cells by Horner.
MUL_PER_SAMPLE = 9
ADD_PER_SAMPLE = 9
# Residual verifier additions the structural terms above do not name individually (index-bit
# subtractions in the sparse selector evaluation vary with the sampled indices). Fitted; the
# resulting model is within 0.25% of measurement on every configuration tested.
ADD_CONST = 930


# =============================================================================
# Security: the QA-code distance bound and the query count
# =============================================================================

def qa_gp(delta, field_bits):
    """The `g_p(delta)` term of the QA distance bound."""
    assert 0.0 < delta < 1.0
    return 1.0 - delta + (delta * math.log2(delta) + (1.0 - delta) * math.log2(1.0 - delta)) / field_bits


def _log2_add(a, b):
    m = max(a, b)
    if math.isinf(m):
        return m
    return m + math.log2(2.0 ** (a - m) + 2.0 ** (b - m))


def qa_distance_failure_log2(delta, row_log_size, c, field_bits):
    """log2 Pr[ delta(C_QA) < delta ] for a random rank-1 systematic QA code of index `c`
    and message length 2**row_log_size. +inf when the bound is vacuous."""
    assert c >= 2
    log_n, log_p = float(row_log_size), float(field_bits)
    eps = qa_gp(delta, field_bits) - (1.0 + log_n / log_p) / c
    if eps <= 0.0:
        return float("inf")

    x = log_p * c * eps
    denom_log = math.log2(-math.expm1(-math.log(2) * x)) if x < 60.0 else 0.0

    t1a = math.log2(c * (c - 1) / 2.0) + log_n - 2.0 * log_p
    t1b = -log_p * math.ceil((c - 1) / (c * delta)) * c * eps - denom_log - log_p
    t2a = math.log2(c) + log_n - log_p
    t2b = -log_p * math.ceil(1.0 / delta) * c * eps - denom_log - log_p
    return min(_log2_add(t1a, t1b), _log2_add(t2a, t2b))


def qa_distance_lower_bound(row_log_size, c, field_bits, failure_bits, iters=100):
    """Largest delta whose failure probability is <= 2**-failure_bits (bisection).
    Returns ~0 when no admissible delta exists (e.g. 31-bit fields)."""
    lo, hi = 1e-9, 1.0 - 1e-9
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if qa_distance_failure_log2(mid, row_log_size, c, field_bits) <= -float(failure_bits):
            lo = mid
        else:
            hi = mid
    return lo


def qa_num_queries(delta, sec_bits):
    """Columns to open: t = ceil(sec_bits / -log2(1 - delta/3)), the provable delta/3 radius."""
    assert 0.0 < delta < 1.0
    return math.ceil(sec_bits / -math.log2(1.0 - delta / 3.0))


def quasar_security(n2, c, field_bits, sec_bits, dist_failure, t, bf_queries):
    """Every term of the paper's Equation (2), in bits."""
    delta = qa_distance_lower_bound(n2, c, field_bits, dist_failure)
    if delta < 1e-6:
        return None
    L = c * (1 << n2)
    lam_dist = -qa_distance_failure_log2(delta, n2, c, field_bits)
    lam_samp = -t * math.log2(1.0 - delta / 3.0)
    lam_aux = -(1.0 + t * math.log2(1.0 - 2.0 * delta / 3.0))
    numerator = 12 * n2 + 2 * c - 4 + (t - 1 + 4 * math.ceil(math.log2(L))) + (delta * L / 3.0 + 1)
    lam_fld = field_bits - math.log2(numerator)
    # BaseFold, rate 1/2 unique decoding: (3/4)^reps. The paper's 249 reps give 103.34 bits.
    lam_bf = bf_queries * -math.log2(0.75)
    lam_vc = 128.0
    total = -math.log2(sum(2.0 ** -b for b in (lam_dist, lam_samp, lam_aux, lam_fld, lam_bf, lam_vc)))
    return dict(delta=delta, t=t, lam_dist=lam_dist, lam_samp=lam_samp, lam_aux=lam_aux,
                lam_fld=lam_fld, lam_bf=lam_bf, lam_total=total)


# =============================================================================
# Transcript: an exact simulation of the duplex sponge
# =============================================================================

class Duplex:
    """DuplexChallenger over the small field: `observe` buffers and permutes every `rate`
    absorbed elements; `sample` permutes whenever the input buffer is dirty or outputs ran out.
    Mirrors `challenger/src/duplex_challenger.rs` exactly, so the permutation count is not a
    fit -- it is the same state machine."""

    def __init__(self, rate):
        self.rate = rate
        self.inp = 0
        self.out = 0
        self.perms = 0

    def observe(self, n=1):
        for _ in range(n):
            self.out = 0
            self.inp += 1
            if self.inp == self.rate:
                self.perms += 1
                self.inp = 0
                self.out = self.rate

    def sample(self, n=1):
        for _ in range(n):
            if self.inp > 0 or self.out == 0:
                self.perms += 1
                self.inp = 0
                self.out = self.rate
            self.out -= 1


def sponge_perms(s, h):
    """Replay `QuasarPcs::verify_core`'s transcript. `s` is the shape dict from `shape()`."""
    limbs, digest = h["limbs"], h["digest_elems"]
    d = Duplex(h["rate"])
    big_obs = lambda n=1: d.observe(limbs * n)   # observe n ed25519 elements
    big_smp = lambda n=1: d.sample(limbs * n)    # sample n ed25519 elements

    def sumcheck(rounds, degree):
        big_obs()                                # claimed_sum
        for _ in range(rounds):
            big_obs(degree)                      # the round polynomial's `degree` evaluations
            big_smp()                            # the round challenge

    def basefold_opening():
        for i in range(s["d"]):
            big_obs(3)                           # h(0), h(1), h(2)
            big_smp()                            # fold challenge
            if i + 1 < s["d"]:
                d.observe(digest)                # folded-oracle root
        big_obs()                                # final value
        d.sample(s["bf_queries"])                # query indices: one small element each

    d.observe(digest)                            # commitment
    d.observe(digest)                            # preprocessed coefficient commitment
    big_smp(s["N1"])                             # proximity weights r
    d.observe(digest)                            # blocks commitment
    d.sample(s["t"])                             # column indices
    big_smp(2)                                   # tau, sigma
    sumcheck(s["log_L"], 2)                      # selector
    big_obs(s["k_on"])                           # block evals at e1
    big_smp(s["n2"] + 2)                         # gamma, alpha, beta
    big_obs(2 * (s["c"] - 1))                    # parity evals at gamma
    sumcheck(s["n2"], 3)                         # encoding, outer layer
    big_obs(s["k_co"])                           # coefficient evals at e2
    sumcheck(s["n2"], 2)                         # encoding, inner layer
    big_obs(2)                                   # input evals at e3
    big_smp()                                    # rho
    sumcheck(s["n2"], 2)                         # claim reduction
    big_obs(s["k_on"] + s["k_co"])               # evals at zeta
    big_smp()                                    # eta
    basefold_opening()                           # online blocks
    basefold_opening()                           # preprocessed coefficients
    return d.perms


# =============================================================================
# Cost model
# =============================================================================

def shape(cfg):
    """Derived shape, plus the security-derived column count."""
    n1, c = cfg["log_rows"], cfg["index"]
    n2 = cfg["num_vars"] - n1
    assert n2 >= 1, "need at least two columns"
    assert c >= 2 and (c & (c - 1)) == 0, "QA index must be a power of two >= 2"
    L = c << n2
    delta = qa_distance_lower_bound(n2, c, cfg["field_bits"], cfg["dist_failure"])
    t = min(qa_num_queries(delta, cfg["sec_bits"]), L) if delta > 1e-6 else None
    n0 = cfg["bf_inv_rate"] << n2                      # BaseFold codeword length
    return dict(
        n1=n1, n2=n2, c=c, N1=1 << n1, N2=1 << n2, L=L,
        log_L=n2 + c.bit_length() - 1,
        t=t, delta=delta,
        d=n2,                                          # BaseFold fold rounds
        n0=n0, log_n0=n0.bit_length() - 1,
        bf_queries=min(cfg["bf_queries"], n0),
        k_on=2 * c,                                    # online blocks committed
        k_co=c - 1,                                    # preprocessed coefficient vectors
    )


def estimate(cfg):
    s = shape(cfg)
    h = HASHES[cfg["hash"]]
    limbs, rate, digest = h["limbs"], h["rate"], h["digest_elems"]
    if s["t"] is None:
        return dict(secure=False, note="no admissible QA distance for this field/length")

    n2, c, t, d = s["n2"], s["c"], s["t"], s["d"]
    qbf = s["bf_queries"]
    leaf_perms = lambda elems: math.ceil(limbs * elems / rate)

    # ---- verifier hashing (exact: one term per hash call in `verify_core`) ----
    # Opened columns: one leaf hash of N1 elements + one path per query.
    col_leaf = t * leaf_perms(s["N1"])
    col_path = t * s["log_L"]
    # A BaseFold opening: per query, a round-0 leaf of 2k values then d-1 leaves of a pair,
    # with paths of shrinking depth.
    path_sum = d * s["log_n0"] - d * (d + 1) // 2
    bf_leaf = lambda k: qbf * (leaf_perms(2 * k) + (d - 1) * leaf_perms(2))
    hash_leaf = col_leaf + bf_leaf(s["k_on"]) + bf_leaf(s["k_co"])
    hash_compress = col_path + 2 * qbf * path_sum
    hash_sponge = sponge_perms(s, h)

    # ---- verifier field arithmetic ----
    # Every term below names a specific piece of `verify_core`, with per-primitive costs
    # measured by the `ed25519_primitives` bin. Nothing here is fitted except ADD_CONST.
    log_c = c.bit_length() - 1
    n1, N1 = s["n1"], s["N1"]
    claims = 5 * c                      # 2c at e1, 2(c-1) at gamma, 2 at e3, 1 at z_lo, c-1 at e2
    # Challenges the verifier draws from the transcript, each recombined by Horner.
    big_samples = N1 + s["log_L"] + 4 * n2 + 2 * d + 6
    # Two inversions per sumcheck round (the barycentric differences and the claim reduction),
    # one per interpolator construction, one per BaseFold round.
    quasar_rounds = s["log_L"] + 3 * n2
    v_inv = 2 * quasar_rounds + 4 + 2 * d
    v_sq = SQ_PER_INV * v_inv

    sc_mul = (INTERP_NEW[2][0] + s["log_L"] * INTERP_EVAL[2][0]      # selector
              + INTERP_NEW[3][0] + n2 * INTERP_EVAL[3][0]            # encoding, outer
              + 2 * (INTERP_NEW[2][0] + n2 * INTERP_EVAL[2][0]))     # inner + reduction
    sc_add = (INTERP_NEW[2][1] + s["log_L"] * INTERP_EVAL[2][1]
              + INTERP_NEW[3][1] + n2 * INTERP_EVAL[3][1]
              + 2 * (INTERP_NEW[2][1] + n2 * INTERP_EVAL[2][1]))
    # Per BaseFold opening: replay the fold sumcheck, then per query recombine the round-0
    # batched pair and walk the fold chain (4 muls per `fold_pair`).
    bf_mul = lambda k: (d * BF_INTERP_MUL + (3 * d - 1) + cfg["bf_inv_rate"]
                        + qbf * (2 * k + 4 * d))
    bf_add = lambda k: d * BF_INTERP_ADD + 3 * d + qbf * (2 * max(k - 1, 0) + 4 * d)

    v_mul = (2 * (N1 - 1) + N1                       # eq(z_hi) table, then r + sigma*z_left
             + t + t * (N1 + 1)                      # tau powers, then one dot per column
             + t * s["log_L"]                        # sparse selector at the terminal point
             + 2 * (c - 1) + 2 * c + 1               # block selector and its combination
             + sc_mul
             + 3 * (c - 1) + 2 * (c - 1)             # batched parity claim, coefficients at e2
             + 2 * (3 * n2 - 1) + 4                  # two H-tilde evaluations and their checks
             + 3 * claims + claims * (3 * n2 + 1)    # claim batching and eq(point, zeta)
             + s["k_on"] + s["k_co"] + (s["k_on"] + s["k_co"])   # eta powers and the two dots
             + bf_mul(s["k_on"]) + bf_mul(s["k_co"])
             + MUL_PER_SAMPLE * big_samples
             + 2)
    v_add = (n1 + N1 + t * N1 + t
             + t * (s["log_L"] // 2 + 1)
             + sc_add
             + log_c + 2 * c + (c - 1) + (c - 2 if c > 2 else 0) + 2 * n2
             + 2 * claims + claims * 3 * n2
             + bf_add(s["k_on"]) + bf_add(s["k_co"])
             + ADD_PER_SAMPLE * big_samples
             + ADD_CONST)

    # ---- proof size ----
    # Field elements: the claimed value; the four sumchecks' claims and round polynomials; the
    # claimed evaluations; the opened columns; and each BaseFold opening's sumcheck, final value
    # and opened rows.
    sumcheck_elems = (1 + 2 * s["log_L"]) + (1 + 3 * n2) + (1 + 2 * n2) + (1 + 2 * n2)
    claim_elems = s["k_on"] + 2 * (c - 1) + s["k_co"] + 2 + (s["k_on"] + s["k_co"])
    bf_elems = lambda k: (3 * d + 1) + qbf * (2 * k + 2 * (d - 1))
    field_elems = (1 + sumcheck_elems + claim_elems + t * s["N1"]
                   + bf_elems(s["k_on"]) + bf_elems(s["k_co"]))
    # Digests: the online blocks commitment, the column paths, and each opening's fold roots and
    # query paths.
    digests = 1 + t * s["log_L"] + 2 * ((d - 1) + qbf * path_sum)
    # Postcard writes a varint length prefix per Vec; one per opened column and per opened row.
    vecs = t * 3 + 2 * qbf * (1 + d * 3) + 20
    proof_bytes = round(field_elems * FIELD_BYTES + digests * h["digest_bytes"] + vecs * VEC_FRAMING)

    # ---- prover operation counts (paper Table 1 + the evaluation phase) ----
    n = 1 << cfg["num_vars"]
    commit_mults = (c - 1) * n
    commit_adds = 2 * (c - 1) * n * n2
    commit_leaf_hashes = s["L"] * leaf_perms(s["N1"])
    commit_compressions = s["L"] - 1
    # Evaluation phase: two row combinations over the whole matrix, two QA encodings of a row,
    # the BaseFold commitments, and the sumchecks (all linear in N2).
    eval_mults = 2 * n + 2 * (c - 1) * s["N2"] + 6 * s["k_on"] * s["N2"]

    sec = quasar_security(n2, c, cfg["field_bits"], cfg["sec_bits"], cfg["dist_failure"], t,
                          qbf if cfg["hash"] != "paper" else cfg["bf_queries"])

    return dict(
        num_vars=cfg["num_vars"], log_rows=s["n1"], log_cols=n2, index=c, N=n,
        delta=round(s["delta"], 6), t=t, bf_queries=qbf,
        sec_bits=round(-t * math.log2(1.0 - s["delta"] / 3.0), 1),
        lam_total=round(sec["lam_total"], 1) if sec else None,
        secure=bool(sec and sec["lam_total"] >= 100.0),
        proof_bytes=proof_bytes, proof_kb=round(proof_bytes / 1000, 1),
        proof_mb=round(proof_bytes / 1e6, 2),
        field_elems=field_elems, digests=digests,
        hash_leaf=hash_leaf, hash_compress=hash_compress, hash_sponge=hash_sponge,
        hash_total=hash_leaf + hash_compress + hash_sponge,
        mul=v_mul, sq=v_sq, inv=v_inv, add=v_add,
        commit_mults=commit_mults, commit_adds=commit_adds,
        commit_hashes=commit_leaf_hashes + commit_compressions,
        eval_mults=eval_mults,
    )


# =============================================================================
# Validation against the profiler
# =============================================================================

# Quantities that are combinatorial and must match exactly, vs. those carrying a fitted
# coefficient and allowed a small relative error.
EXACT_KEYS = ["t", "hash_leaf", "hash_compress", "hash_sponge", "inv", "sq", "mul"]
APPROX_KEYS = {"add": 1.0, "proof_bytes": 2.0}


def check(measured_path):
    rows = json.load(open(measured_path))
    print(f"{'config':>12} {'quantity':>14} {'predicted':>12} {'measured':>12} {'err':>8}")
    worst = {}
    failures = 0
    for row in rows:
        cfg = dict(CONFIG)
        cfg.update(num_vars=row["num_vars"], log_rows=row["log_rows"], index=row["index"],
                   bf_inv_rate=row["bf_inv_rate"], bf_queries=row["bf_queries"])
        est = estimate(cfg)
        name = f"{row['num_vars']}:{row['log_rows']}:{row['index']}"
        for key in EXACT_KEYS + list(APPROX_KEYS):
            p, m = est[key], row["q" if key == "t" else key]
            err = 0.0 if m == p else (100.0 * (p - m) / m if m else float("inf"))
            ok = (p == m) if key in EXACT_KEYS else abs(err) <= APPROX_KEYS[key]
            worst[key] = max(worst.get(key, 0.0), abs(err))
            if not ok:
                failures += 1
                print(f"{name:>12} {key:>14} {p:>12} {m:>12} {err:>7.2f}%  MISMATCH")
    print()
    print(f"{'quantity':>14} {'worst |err|':>12}   verdict")
    for key in EXACT_KEYS:
        print(f"{key:>14} {worst[key]:>11.2f}%   {'EXACT' if worst[key] == 0 else 'MISMATCH'}")
    for key, tol in APPROX_KEYS.items():
        print(f"{key:>14} {worst[key]:>11.2f}%   {'ok' if worst[key] <= tol else 'OVER'} (tol {tol}%)")
    print()
    print(f"{len(rows)} configurations, {failures} mismatches")
    return failures == 0


# ---- paper Table 4: (k, n2, delta_QA, t) at Mersenne127, c = 2, lambda_dist = 103 ----
TABLE_4 = [
    (20, 12, 0.444951, 445), (21, 13, 0.441024, 450), (22, 14, 0.437098, 454),
    (23, 15, 0.433173, 458), (24, 16, 0.429248, 463), (25, 17, 0.425324, 467),
    (26, 18, 0.421401, 472), (27, 19, 0.417478, 477), (28, 20, 0.413556, 482),
    (29, 21, 0.409635, 487),
]


def self_test():
    print("paper Table 4 (Mersenne127, c=2, lambda_dist=103, 249 BaseFold reps)")
    print(f"{'k':>3} {'n2':>3} {'delta_paper':>12} {'delta_here':>12} {'t_paper':>8} {'t_here':>7} "
          f"{'lam_total':>10}  ok")
    ok_all = True
    for k, n2, d_paper, t_paper in TABLE_4:
        d = qa_distance_lower_bound(n2, 2, 127, 103)
        t = qa_num_queries(d, 103)
        sec = quasar_security(n2, 2, 127, 103, 103, t, 249)
        ok = (round(d, 6) == d_paper) and (t == t_paper) and abs(sec["lam_total"] - 101.5) < 0.15
        ok_all &= ok
        print(f"{k:>3} {n2:>3} {d_paper:>12.6f} {d:>12.6f} {t_paper:>8} {t:>7} "
              f"{sec['lam_total']:>10.1f}  {'OK' if ok else 'MISMATCH'}")
    print("Table 4:", "ALL ROWS MATCH" if ok_all else "MISMATCH")
    return ok_all


def estimate_stark(log_height, widths, overrides=None):
    """Whole-proof estimate for a multi-stark AIR: one Quasar instance per committed trace,
    each sized for that trace's stacked arity `log_height + ceil(log2 width)`.

    The multi-stark framework's own contribution (zerocheck transcripts, per-column opening
    evaluations, the layout sumcheck) is O(width + log_height) field elements — under 1% of a
    PCS-dominated proof, so it is NOT modelled. Measured against the ed25519 examples:

        fib arity 11 / 15 / 18   +0.9% / +0.7% / +0.7%
        eddsa (2 instances)      -2.3%
        fib arity 7 (15 KB)      -7.8%   <- framework overhead visible only when the proof is tiny

    `widths` is one entry per committed trace (main first, then preprocessed if the AIR has one).
    """
    cfg = dict(CONFIG)
    cfg.update(overrides or {})
    instances = []
    for width in widths:
        stacked = log_height + max(1, math.ceil(math.log2(width)))
        shape = choose_shape(stacked, cfg)
        if shape is None:
            return None
        per = dict(cfg)
        per.update(num_vars=stacked, log_rows=shape)
        instances.append(estimate(per))
    return dict(
        log_height=log_height,
        widths=list(widths),
        instances=instances,
        proof_bytes=sum(i["proof_bytes"] for i in instances),
        hash_total=sum(i["hash_total"] for i in instances),
        mul=sum(i["mul"] for i in instances),
        secure=all(i["secure"] for i in instances),
    )


def choose_shape(num_vars, cfg):
    """The `log_rows` minimising predicted proof size — the Python twin of
    `QuasarParams::choose` / `choose_best_effort`, pinned to it by
    `chooser_matches_cost_model` on the Rust side.

    Prefers shapes whose codeword can hold the unclamped query count; if none can,
    falls back to the smallest-proof shape at whatever soundness the clamp allows.
    """
    best = best_effort = None
    for log_rows in range(num_vars):
        n2 = num_vars - log_rows
        delta = qa_distance_lower_bound(n2, cfg["index"], cfg["field_bits"], cfg["dist_failure"])
        if delta < 1e-6:
            continue
        reaches_target = qa_num_queries(delta, cfg["sec_bits"]) <= cfg["index"] << n2
        per = dict(cfg)
        per.update(num_vars=num_vars, log_rows=log_rows)
        bytes_ = estimate(per)["proof_bytes"]
        if best_effort is None or bytes_ < best_effort[1]:
            best_effort = (log_rows, bytes_)
        if reaches_target and (best is None or bytes_ < best[1]):
            best = (log_rows, bytes_)
    chosen = best or best_effort
    return chosen[0] if chosen else None


def report(est):
    if not est.get("secure", True) and "note" in est:
        print(est["note"])
        return
    print(f"N = 2^{est['num_vars']} = 2^{est['log_rows']} rows x 2^{est['log_cols']} cols, "
          f"c = {est['index']}, delta_QA = {est['delta']:.6f}")
    print(f"  security   : {est['t']} columns opened -> {est['sec_bits']} bits sampling, "
          f"{est['lam_total']} bits total ({'>=100' if est['secure'] else 'BELOW 100'})")
    print(f"  proof      : {est['proof_mb']} MB  ({est['field_elems']} field elements, "
          f"{est['digests']} digests)")
    print(f"  verifier   : {est['hash_total']} perms "
          f"(leaf {est['hash_leaf']}, compress {est['hash_compress']}, sponge {est['hash_sponge']})")
    print(f"               {est['mul']} mul, {est['sq']} sq, {est['inv']} inv, {est['add']} add")
    print(f"  commit     : {est['commit_mults']} mults, {est['commit_adds']} adds, "
          f"{est['commit_hashes']} hash perms")
    print(f"  eval prover: ~{est['eval_mults']} mults")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("params", nargs="?", help="JSON file of parameter overrides (dict or list)")
    ap.add_argument("--json", action="store_true", help="machine-readable JSON output")
    ap.add_argument("--check", metavar="MEASURED", help="diff predictions against profiler output")
    ap.add_argument(
        "--stark",
        metavar="LOG_HEIGHT:WIDTH[,WIDTH...]",
        help="whole-proof estimate for a multi-stark AIR, e.g. --stark 8:45,8 for the EdDSA AIR",
    )
    args = ap.parse_args()

    if args.check:
        sys.exit(0 if check(args.check) else 1)

    if args.stark:
        log_height, widths = args.stark.split(":")
        widths = [int(w) for w in widths.split(",")]
        est = estimate_stark(int(log_height), widths)
        if est is None:
            print("no admissible shape")
            sys.exit(1)
        if args.json:
            print(json.dumps(est, indent=2))
        else:
            print(f"trace 2^{est['log_height']} rows, committed widths {est['widths']}")
            for width, inst in zip(est["widths"], est["instances"]):
                print(f"  width {width:<4} -> 2^{inst['num_vars']} stacked, "
                      f"N1=2^{inst['log_rows']}, t={inst['t']}, {inst['proof_mb']} MB, "
                      f"{inst['lam_total']} bits")
            print(f"  total proof {est['proof_bytes'] / 1e6:.2f} MB, "
                  f"{est['hash_total']} verifier perms"
                  f"{'' if est['secure'] else '   (BELOW 100 bits)'}")
        sys.exit(0)

    if args.params:
        loaded = json.load(open(args.params))
        cfgs = loaded if isinstance(loaded, list) else [loaded]
    else:
        if not self_test():
            sys.exit(1)
        print()
        cfgs = [{}]

    results = []
    for override in cfgs:
        cfg = dict(CONFIG)
        cfg.update(override)
        results.append(estimate(cfg))

    if args.json:
        print(json.dumps(results if len(results) > 1 else results[0], indent=2))
    else:
        for est in results:
            report(est)
            print()


if __name__ == "__main__":
    main()
