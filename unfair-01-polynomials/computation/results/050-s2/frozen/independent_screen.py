#!/usr/bin/env python3
"""Independent exact necessary mother screen for degree-7 S={2}.

This is intentionally a fresh implementation.  Elements of Q(a) are stored
as U+V*a, with a^2=T*a-N; no floating-point arithmetic is used.  For each
deduplicated (T,N,Delta,b,c) low prefix, the raw ledger is the full lexicographic
box d=-52..52, e=-29..29, f=-8..8, one byte per point.

It applies only necessary tests on P and H=P*Pbar.  A surviving row is not a
factorization, nor a counterexample to the uniform-sum problem.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path


# Seven exact positive arguments p/q, in the order specified for this screen.
POINTS = ((1, 1), (1, 2), (2, 1), (1, 4), (4, 1), (3, 4), (4, 3))
DEG = 7
RAW_PER_CORE = 105 * 59 * 17


def qadd(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (x[0] + y[0], x[1] + y[1])


def qmul(x: tuple[int, int], y: tuple[int, int], T: int, N: int) -> tuple[int, int]:
    """(u+v*a)(r+s*a), with a^2=T*a-N."""
    u, v = x
    r, s = y
    return (u * r - N * v * s, u * s + v * r + T * v * s)


def qbar(x: tuple[int, int], T: int) -> tuple[int, int]:
    """The conjugate of u+v*a, still expressed in the a-basis."""
    u, v = x
    return (u + v * T, -v)


def abs_both_lt(x: tuple[int, int], T: int, D: int, num: int, den: int) -> bool:
    """Strictly |sigma(x)| < num/den for both real embeddings."""
    u, v = x
    # 2 sigma(x) = R +/- v sqrt(D).  Its largest absolute value is
    # |R|+|v|sqrt(D), so this one squared comparison checks both embeddings.
    R = 2 * u + v * T
    room = 2 * num - den * abs(R)
    return room > 0 and (den * v) * (den * v) * D < room * room


def positive_both(x: tuple[int, int], T: int, D: int) -> bool:
    """Strictly sigma(x)>0 in both embeddings."""
    u, v = x
    R = 2 * u + v * T
    return R > 0 and v * v * D < R * R


def newton_sums(p: list[tuple[int, int]], T: int, N: int) -> list[tuple[int, int]]:
    """Derive s_1,...,s_7 from k p_k+sum_{i=1}^k p_i s_{k-i}=0."""
    s = [(DEG, 0)]  # s_0 is not used below; position makes indexing explicit.
    for k in range(1, DEG + 1):
        total = (k * p[k][0], k * p[k][1])
        for i in range(1, k):
            total = qadd(total, qmul(p[i], s[k - i], T, N))
        s.append((-total[0], -total[1]))
    return s


def eval_scaled(p: list[tuple[int, int]], num: int, den: int) -> tuple[int, int]:
    """Return den^7 P(num/den), represented in Z[a]."""
    out = (0, 0)
    for i, coeff in enumerate(p):
        scale = num ** i * den ** (DEG - i)
        out = qadd(out, (scale * coeff[0], scale * coeff[1]))
    return out


def norm_poly(p: list[tuple[int, int]], T: int, N: int) -> list[int]:
    """Compute P*Pbar directly in Z[x], asserting exact cancellation."""
    h = [(0, 0) for _ in range(2 * DEG + 1)]
    for i, x in enumerate(p):
        for j, y in enumerate(p):
            h[i + j] = qadd(h[i + j], qmul(x, qbar(y, T), T, N))
    if any(v for _, v in h):
        raise AssertionError("P*Pbar was not rational")
    return [u for u, _ in h]


def disk_failure(ascending: list[int]) -> int | None:
    """First strict-Schur failure degree for roots in the open unit disk."""
    a = list(ascending)
    while len(a) > 1:
        divisor = math.gcd(*a)
        if a[-1] < 0:
            divisor = -divisor
        a = [x // divisor for x in a]
        if abs(a[0]) >= a[-1]:
            return len(a) - 1
        first, last = a[0], a[-1]
        n = len(a) - 1
        a = [last * a[i] - first * a[n - i] for i in range(1, n + 1)]
    return None


def check_candidate(T: int, N: int, D: int, b: int, c: int, d: int, e: int, f: int) -> int:
    """Return the first failure code, or 255.  The caller records all counts."""
    p = [(1, 0), (b, 0), (0, 1), (c, 0), (d, 0), (e, 0), (f, 0), (1, 0)]
    sums = newton_sums(p, T, N)
    for k in range(4, 8):
        # The prefix input has already independently enforced s1..s3.
        uk_num = 3 * (2 ** (2 * k) + 1) + 2 ** k
        uk_den = 2 ** k  # 3(2^k+2^-k)+1
        if not abs_both_lt(sums[k], T, D, uk_num, uk_den):
            return k
    for j, (num, den) in enumerate(POINTS):
        if not positive_both(eval_scaled(p, num, den), T, D):
            return 10 + j
    h = norm_poly(p, T, N)
    # H(2x) has roots lambda/2, hence certifies |lambda|<2.
    failure = disk_failure([coef << i for i, coef in enumerate(h)])
    if failure is not None:
        return 20 + failure
    # (x^14 H(1/x))(2x) has roots 1/(2 lambda), certifying |lambda|>1/2.
    failure = disk_failure([coef << i for i, coef in enumerate(h[::-1])])
    if failure is not None:
        return 40 + failure
    return 255


def atomic_bytes(path: Path, data: bytes) -> None:
    temp = path.with_suffix(path.suffix + ".partial")
    with temp.open("wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temp, path)


def canonical_cores(prefix: dict) -> list[tuple[int, int, int, int, int]]:
    rows = {(r["T"], r["N"], r["Delta"], r["b"], r["c"]) for r in prefix["survivors"]}
    out = sorted(rows)
    if len(out) != 359:
        raise AssertionError(f"expected 359 distinct cores, got {len(out)}")
    return out


def screen_one(core: tuple[int, int, int, int, int]) -> tuple[dict, bytes, list[list[int]]]:
    T, N, D, b, c = core
    counts: Counter = Counter(raw=0)
    ledger = bytearray()
    survivors: list[list[int]] = []
    for d in range(-52, 53):
        for e in range(-29, 30):
            for f in range(-8, 9):
                code = check_candidate(T, N, D, b, c, d, e, f)
                counts["raw"] += 1
                if code == 4:
                    counts["reject4"] += 1
                elif code == 5:
                    counts["pass4"] += 1
                    counts["reject5"] += 1
                elif code == 6:
                    counts["pass4"] += 1; counts["pass5"] += 1
                    counts["reject6"] += 1
                elif code == 7:
                    counts["pass4"] += 1; counts["pass5"] += 1; counts["pass6"] += 1
                    counts["reject7"] += 1
                else:
                    counts["pass4"] += 1; counts["pass5"] += 1; counts["pass6"] += 1; counts["pass7"] += 1
                    if 10 <= code <= 16:
                        counts[f"axis_{code - 10}"] += 1
                    else:
                        counts["positive_axis"] += 1
                        if 20 <= code < 40:
                            counts["outer_rejected"] += 1
                        elif 40 <= code < 60:
                            counts["inner_rejected"] += 1
                        elif code == 255:
                            counts["annulus_survived"] += 1
                            survivors.append([T, N, D, b, c, d, e, f])
                        else:
                            raise AssertionError(f"unknown failure code {code}")
                ledger.append(code)
    if len(ledger) != RAW_PER_CORE:
        raise AssertionError((len(ledger), RAW_PER_CORE))
    return dict(counts), bytes(ledger), survivors


def main() -> None:
    if not __debug__:
        raise RuntimeError("do not run with python -O")
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=359)
    ap.add_argument("--seconds", type=float, default=14400)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    prefix_bytes = args.prefix.read_bytes()
    prefix = json.loads(prefix_bytes)
    cores = canonical_cores(prefix)
    if args.start < 0 or args.count < 0 or args.start + args.count > len(cores):
        raise ValueError("invalid core slice")
    sha = lambda data: hashlib.sha256(data).hexdigest()
    binding = {"source_sha256": sha(Path(__file__).read_bytes()), "prefix_sha256": sha(prefix_bytes)}
    started = time.monotonic()
    done: list[int] = []
    for i in range(args.start, args.start + args.count):
        json_path = args.out / f"core-{i:03d}.json"
        code_path = args.out / f"core-{i:03d}-codes.bin.gz"
        if json_path.exists() and code_path.exists():
            old = json.loads(json_path.read_text())
            if old.get("status") == "complete" and old.get("binding") == binding and old.get("core") == list(cores[i]) and old.get("gzip_code_sha256") == sha(code_path.read_bytes()):
                done.append(i)
                continue
            raise RuntimeError(f"mismatched existing checkpoint {json_path}")
        if time.monotonic() - started > args.seconds:
            break
        one_started = time.monotonic()
        counts, codes, survivors = screen_one(cores[i])
        packed = gzip.compress(codes, mtime=0)
        atomic_bytes(code_path, packed)
        receipt = {
            "status": "complete", "index": i, "core": list(cores[i]), "binding": binding,
            "counts": counts, "raw_code_sha256": sha(codes), "gzip_code_sha256": sha(packed),
            "code_length": len(codes), "survivors": survivors,
            "elapsed_seconds": time.monotonic() - one_started,
        }
        atomic_bytes(json_path, (json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n").encode())
        done.append(i)
        print(json.dumps({"completed": len(done), "index": i, "seconds": receipt["elapsed_seconds"], "counts": counts}, sort_keys=True), flush=True)
    aggregate: Counter = Counter()
    survivors: list[list[int]] = []
    for i in done:
        item = json.loads((args.out / f"core-{i:03d}.json").read_text())
        aggregate.update(item["counts"])
        survivors.extend(item["survivors"])
    status = "complete" if len(done) == args.count else "partial_unknown"
    summary = {
        "status": status, "scope": "Necessary P/H filters only; any annulus survivor needs separate cofactor continuation.",
        "requested_start": args.start, "requested_count": args.count, "completed_indices": done,
        "binding": binding, "raw_per_core": RAW_PER_CORE, "counts": dict(aggregate), "survivors": survivors,
        "elapsed_seconds": time.monotonic() - started, "python": sys.version,
        "code_schema": "Full lex d,e,f box: 4..7=first Newton failure; 10..16=first positive-axis failure at (1,1/2,2,1/4,4,3/4,4/3); 20+degree=outer Schur; 40+degree=inner Schur; 255=survivor.",
    }
    atomic_bytes(args.out / "summary.json", (json.dumps(summary, sort_keys=True, indent=2) + "\n").encode())
    print(json.dumps({"status": status, "completed": len(done), "counts": dict(aggregate), "seconds": summary["elapsed_seconds"]}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
