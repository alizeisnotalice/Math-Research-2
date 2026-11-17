#!/usr/bin/env python3
"""Exact unrestricted-integer-cofactor prefix tree for degree-seven S={2}.

Rows are [T,N,D,b,c,d,e,f].  Arithmetic is exact in Z[a]/(a^2-T*a+N).
The tree deliberately permits signed, arbitrary-degree integer G,Q.  It adds
only the already-proved S={2} necessities G3=Q3=0 and A5=B5=0.  A row is
excluded only if every prefix dies; resource limits remain explicitly unknown.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path


Pair = tuple[int, int]


def add(x: Pair, y: Pair) -> Pair:
    return (x[0] + y[0], x[1] + y[1])


def product(x: Pair, y: Pair, t: int, n: int) -> Pair:
    """(p+qa)(r+sa), where a^2=T*a-N."""
    return (x[0] * y[0] - n * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0] + t * x[1] * y[1])


def floor_pair(x: Pair, t: int, disc: int) -> int:
    """Exact floor of p+q*(T+sqrt(D))/2; D is positive nonsquare."""
    p, q = x
    if q == 0:
        return p
    r = math.isqrt(q * q * disc)
    # For q<0, floor(q sqrt(D))=-floor(|q|sqrt(D))-1, strictly because D
    # is nonsquare.  Dividing after this replacement is exact.
    return (2 * p + q * t + (r if q > 0 else -r - 1)) // 2


def integer_choices(known: Pair, t: int, disc: int) -> tuple[int, ...]:
    """All z in Z for which known+z lies in [0,1]."""
    floor = floor_pair(known, t, disc)
    return (-floor, 1 - floor) if known[1] == 0 else (-floor,)


def before_free(p: tuple[Pair, ...], coeffs: tuple[int, ...], k: int) -> Pair:
    """Coefficient [x^k](P*G) before adding free integer G_k."""
    u = v = 0
    for j in range(1, min(7, k) + 1):
        pj = p[j]
        gj = coeffs[k - j]
        u += pj[0] * gj
        v += pj[1] * gj
    return (u, v)


def convolution_coefficient(a: tuple[Pair, ...], b: tuple[Pair, ...], k: int,
                            t: int, n: int) -> Pair:
    out = (0, 0)
    for i in range(k + 1):
        j = k - i
        if i < len(a) and j < len(b):
            out = add(out, product(a[i], b[j], t, n))
    return out


def digest_states(states: list[tuple]) -> str:
    h = hashlib.sha256()
    for g, q, _, _ in sorted(states, key=lambda z: (z[0], z[1])):
        h.update(json.dumps([g, q], separators=(",", ":")).encode() + b"\n")
    return h.hexdigest()


def validate_row(row: tuple[int, ...]) -> None:
    if len(row) != 8:
        raise ValueError(f"S2 row must have 8 integers, got {row}")
    t, n, disc, *_ = row
    if disc <= 0 or math.isqrt(disc) ** 2 == disc or t * t - 4 * n != disc:
        raise ValueError(f"invalid real-quadratic row {row}")


def run_one(row: tuple[int, ...], depth: int, state_cap: int, seconds: float) -> dict:
    validate_row(row)
    t, n, disc, b, c, d, e, f = row
    p = ((1, 0), (b, 0), (0, 1), (c, 0), (d, 0), (e, 0), (f, 0), (1, 0))
    pb = ((1, 0), (b, 0), (t, -1), (c, 0), (d, 0), (e, 0), (f, 0), (1, 0))
    states: list[tuple] = [((1,), (1,), ((1, 0),), ((1, 0),))]
    trace = []
    start = time.monotonic()
    status = "unknown_depth"

    for k in range(1, depth + 1):
        nxt = []
        rejected = {"forced_endpoint": 0, "forced_g3": 0, "nonbinary_C": 0}
        for g, q, aa, bb in states:
            left = before_free(p, g, k)
            right = before_free(pb, q, k)
            gchoices = integer_choices(left, t, disc)
            qchoices = integer_choices(right, t, disc)
            for gi in gchoices:
                aval = (left[0] + gi, left[1])
                for qi in qchoices:
                    bval = (right[0] + qi, right[1])
                    if k == 3 and (gi != 0 or qi != 0):
                        rejected["forced_g3"] += 1
                        continue
                    if k == 5 and (aval != (0, 0) or bval != (0, 0)):
                        rejected["forced_endpoint"] += 1
                        continue
                    cc = convolution_coefficient(aa + (aval,), bb + (bval,), k, t, n)
                    if cc not in ((0, 0), (1, 0)):
                        rejected["nonbinary_C"] += 1
                        continue
                    nxt.append((g + (gi,), q + (qi,), aa + (aval,), bb + (bval,)))
        unique = {(g, q): (g, q, aa, bb) for g, q, aa, bb in nxt}
        if len(unique) != len(nxt):
            raise AssertionError("duplicate G,Q prefix")
        states = list(unique.values())
        trace.append({"k": k, "live": len(states), "state_sha256": digest_states(states),
                      "rejected": rejected})
        if not states:
            status = "dead_all_prefixes"
            break
        if len(states) > state_cap:
            status = "unknown_state_cap"
            break
        if time.monotonic() - start > seconds:
            status = "unknown_time_limit"
            break

    out = {"row": list(row), "status": status, "depth_reached": len(trace), "trace": trace,
           "elapsed_seconds": time.monotonic() - start,
           "limits": {"depth": depth, "state_cap": state_cap, "seconds": seconds},
           "constraints": {"G3_Q3_zero": True, "A5_B5_zero": True}}
    if status != "dead_all_prefixes":
        out["live_prefixes"] = [{"G": list(g), "Q": list(q)} for g, q, _, _ in states]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path, help="JSON rows [T,N,D,b,c,d,e,f] or object with rows")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int)
    ap.add_argument("--depth", type=int, default=48)
    ap.add_argument("--state-cap", type=int, default=25000)
    ap.add_argument("--seconds-per-core", type=float, default=8.0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    raw = args.input.read_bytes()
    parsed = json.loads(raw)
    rows = parsed["rows"] if isinstance(parsed, dict) else parsed
    stop = len(rows) if args.stop is None else min(args.stop, len(rows))
    records = []
    for index, row in enumerate(rows[args.start:stop], args.start):
        rec = run_one(tuple(row), args.depth, args.state_cap, args.seconds_per_core)
        rec["input_index"] = index
        records.append(rec)
        if not args.quiet:
            print(json.dumps({"index": index, "status": rec["status"], "depth": rec["depth_reached"],
                              "seconds": rec["elapsed_seconds"]}), flush=True)
    counts: dict[str, int] = {}
    for rec in records:
        counts[rec["status"]] = counts.get(rec["status"], 0) + 1
    args.out.write_text(json.dumps({"input_sha256": hashlib.sha256(raw).hexdigest(),
                                    "range": [args.start, stop], "records": records,
                                    "status_counts": counts}, indent=2) + "\n")


if __name__ == "__main__":
    main()
