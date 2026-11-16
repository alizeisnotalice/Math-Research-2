#!/usr/bin/env python3
"""Independent audit for r022 terminal tau=2 lift traces.

It deliberately recomputes root pairs, earliest low-root-neighbour rank,
L, L*d^2, L*F0, L*G0, group sums, and the exact quadratic identity from the
saved terminal prefixes; it does not import the Swift/Metal implementation.
When an old same-seed JSONL is supplied, all pre-existing result fields are
also compared (timings excluded).
"""
import json
import sys


def edge(a, b, n):
    return abs(a % n - b % n) + abs(a // n - b // n) == 1


def ints(values):
    return [int(x) for x in values]


def audit(new_jsonl):
    rows = [json.loads(x) for x in open(new_jsonl, encoding="utf-8")]
    meta = rows[0]
    rep = next(x for x in rows if x.get("kind") == "replicate")
    assert meta["n"] * meta["n"] == len(json.load(open(new_jsonl + ".trace.json"))["rows"])
    trace = json.load(open(new_jsonl + ".trace.json", encoding="utf-8"))
    assert trace["tau"] == 2 and trace["P"] == meta["particles"]
    term = trace["rows"][-1]
    n, N, P = meta["n"], meta["n"] ** 2, meta["particles"]
    assert term["terminal"] and term["rank"] == N
    groups = []
    for g in range((P + 255) // 256):
        acc = [0, 0, 0, 0]
        for pref, state in zip(term["prefix"][256 * g:min(P, 256 * (g + 1))], term["state"][256 * g:min(P, 256 * (g + 1))]):
            assert len(pref) == N - 1
            x, y = pref[:2]
            assert (state[2], state[3]) == (x, y)
            assert 0 <= x < N and 0 <= y < N and x != y and not edge(x, y, n)
            rank = {v: k + 1 for k, v in enumerate(pref)}
            assert len(rank) == N - 1 and y not in [v for v in pref[2:] if v == y]
            j = min(rank.get(v, N) for v in range(N) if edge(y, v, n))
            assert 3 <= j <= N - 1
            L = j - 2
            d2 = (x % n - y % n) ** 2 + (x // n - y // n) ** 2
            ax, ay = 2 * (x % n) - (n - 1), 2 * (x // n) - (n - 1)
            bx, by = 2 * (y % n) - (n - 1), 2 * (y // n) - (n - 1)
            f0 = ax * ax + ay * ay + bx * bx + by * by
            g0 = ax * bx + ay * by
            assert L * f0 - 2 * L * g0 == 4 * L * d2
            acc[0] += L; acc[1] += L * d2; acc[2] += L * f0; acc[3] += L * g0
        groups.append(acc)
    keys = ["tau2_lift_L_groups", "tau2_lift_L_distance2_groups", "tau2_lift_L_F0_groups", "tau2_lift_L_G0_groups"]
    observed = [ints(rep[k]) for k in keys]
    expected = [[x[i] for x in groups] for i in range(4)]
    assert observed == expected, (observed, expected)
    residual = [x[2] - 2 * x[3] - 4 * x[1] for x in groups]
    assert ints(rep["tau2_lift_group_LF_minus_2LG_minus_4LD"]) == residual == [0] * len(groups)
    assert term["tau2_lift_bad"] == 0
    return {"n": n, "P": P, "groups": groups, "raw": rep}


def same_seed_old(new_jsonl, old_jsonl):
    new = next(x for x in map(json.loads, open(new_jsonl)) if x.get("kind") == "replicate")
    old = next(x for x in map(json.loads, open(old_jsonl)) if x.get("kind") == "replicate")
    excluded = {"gpu_seconds", "wall_seconds"}
    common = (set(new) & set(old)) - excluded
    mismatch = {k: (new[k], old[k]) for k in common if new[k] != old[k]}
    assert not mismatch, mismatch
    return sorted(common)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(f"usage: {sys.argv[0]} NEW.jsonl [OLD.jsonl]")
    got = audit(sys.argv[1])
    common = same_seed_old(sys.argv[1], sys.argv[2]) if len(sys.argv) == 3 else []
    print(json.dumps({"ok": True, "n": got["n"], "P": got["P"], "groups": got["groups"], "old_common_fields_identical": common}, sort_keys=True))
