#!/usr/bin/env python3
"""Exact unrestricted-integer-cofactor prefix tree for degree-seven S={3} cores.

Input rows have [T,N,D,c,d,e,f].  All arithmetic is in Z[a]/(a^2-T*a+N).
A core is excluded only when its entire tree is empty.  Depth, state, or time
limits are reported as unknown rather than exclusions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path


def add(x, y): return (x[0] + y[0], x[1] + y[1])


def floor_pair(x, t, disc):
    p, q = x
    if q == 0:
        return p
    rad = math.isqrt(q * q * disc)
    return (2 * p + q * t + (rad if q > 0 else -rad - 1)) // 2


def choices(base, t, disc):
    p, q = base
    z = floor_pair(base, t, disc)
    return (-z, 1 - z) if q == 0 else (-z,)


def base(p, coeffs, k):
    # Coefficient of P*G at k before the free integer G_k.
    r = s = 0
    for j in range(1, min(k, 7) + 1):
        u, v = p[j]
        z = coeffs[k - j]
        r += u * z
        s += v * z
    return (r, s)


def ccoef(a, b, k):
    r = s = 0
    for i in range(k + 1):
        u, v = a[i]
        x, y = b[k - i]
        # Product in a^2=T*a-N is never needed here: C must be rational,
        # but multiplication parameters are supplied by closure below.
        r += u * x
        s += u * y + v * x
    return r, s


def product(x, y, t, n):
    return (x[0] * y[0] - n * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0] + t * x[1] * y[1])


def ccoef_field(a, b, k, t, n):
    out = (0, 0)
    for i in range(k + 1):
        out = add(out, product(a[i], b[k - i], t, n))
    return out


def digest_states(states):
    h = hashlib.sha256()
    for g, q, _, _ in sorted(states, key=lambda z: (z[0], z[1])):
        h.update(json.dumps([g, q], separators=(",", ":")).encode() + b"\n")
    return h.hexdigest()


def run_one(row, depth, state_cap, seconds):
    t, n, disc, c, d, e, f = row
    p = ((1, 0), (0, 0), (c, 0), (0, 1), (d, 0), (e, 0), (f, 0), (1, 0))
    pb = ((1, 0), (0, 0), (c, 0), (t, -1), (d, 0), (e, 0), (f, 0), (1, 0))
    states = [((1,), (1,), ((1, 0),), ((1, 0),))]
    trace = []
    start = time.monotonic()
    status = "unknown_depth"
    for k in range(1, depth + 1):
        nxt = []
        rejected = {"forced_endpoint": 0, "nonbinary_C": 0}
        for g, q, aa, bb in states:
            left, right = base(p, g, k), base(pb, q, k)
            for gi in choices(left, t, disc):
                aval = (left[0] + gi, left[1])
                for qi in choices(right, t, disc):
                    bval = (right[0] + qi, right[1])
                    if k in (1, 4) and (aval != (0, 0) or bval != (0, 0)):
                        rejected["forced_endpoint"] += 1
                        continue
                    cc = ccoef_field(aa + (aval,), bb + (bval,), k, t, n)
                    if cc != (0, 0) and cc != (1, 0):
                        rejected["nonbinary_C"] += 1
                        continue
                    nxt.append((g + (gi,), q + (qi,), aa + (aval,), bb + (bval,)))
        unique = {(g, q): (g, q, aa, bb) for g, q, aa, bb in nxt}
        if len(unique) != len(nxt):
            raise AssertionError("duplicate G,Q prefix")
        states = list(unique.values())
        trace.append({"k": k, "live": len(states), "state_sha256": digest_states(states), "rejected": rejected})
        if not states:
            status = "dead_all_prefixes"
            break
        if len(states) > state_cap:
            status = "unknown_state_cap"
            break
        if time.monotonic() - start > seconds:
            status = "unknown_time_limit"
            break
    result = {
        "row": list(row), "status": status, "depth_reached": len(trace),
        "trace": trace, "elapsed_seconds": time.monotonic() - start,
        "limits": {"depth": depth, "state_cap": state_cap, "seconds": seconds},
    }
    if status != "dead_all_prefixes":
        result["live_prefixes"] = [{"G": list(g), "Q": list(q)} for g, q, _, _ in states]
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path, help="JSON list of [T,N,D,c,d,e,f]")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int)
    ap.add_argument("--depth", type=int, default=48)
    ap.add_argument("--state-cap", type=int, default=25000)
    ap.add_argument("--seconds-per-core", type=float, default=8.0)
    args = ap.parse_args()
    raw_input = args.input.read_bytes()
    parsed = json.loads(raw_input)
    rows = parsed["rows"] if isinstance(parsed, dict) else parsed
    stop = len(rows) if args.stop is None else min(args.stop, len(rows))
    selected = rows[args.start:stop]
    records = []
    for i, row in enumerate(selected, args.start):
        rec = run_one(tuple(row), args.depth, args.state_cap, args.seconds_per_core)
        rec["input_index"] = i
        records.append(rec)
        print(json.dumps({"index": i, "status": rec["status"], "depth": rec["depth_reached"], "seconds": rec["elapsed_seconds"]}), flush=True)
    summary = {"input_sha256": hashlib.sha256(raw_input).hexdigest(),
               "range": [args.start, stop], "records": records,
               "status_counts": {}}
    for rec in records:
        summary["status_counts"][rec["status"]] = summary["status_counts"].get(rec["status"], 0) + 1
    args.out.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
