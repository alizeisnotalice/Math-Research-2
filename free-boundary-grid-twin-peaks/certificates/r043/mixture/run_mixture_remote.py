#!/usr/bin/env python3
"""Finite, reproducible PRNG diagnostic for the exact curved-event sampler.

The curve draw is performed here with Python's integer ``randrange(S1)``;
the C++ binary only samples exactly conditional on a fixed curve.  Every
emitted word is checked both for its nominated event and against every binary
curve when computing X.  Thus a failed C++ job, malformed word, or structural
failure makes the whole run fail closed.

This program produces *finite numerical evidence*.  It neither estimates the
original conditional law directly nor certifies an asymptotic statement.
"""
import argparse
import concurrent.futures as futures
import hashlib
import json
import os
import random
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

SCHEMA = "curved-mixture-prng-diagnostic-v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def adjacency(n):
    return [sum(1 << u for u in range(n * n)
                if abs(u % n - v % n) + abs(u // n - v // n) == 1)
            for v in range(n * n)]


def curve_data(n, h, adj):
    H = sum(1 << (y * n + x) for y, c in enumerate(h) for x in range(int(c)))
    full = (1 << (n * n)) - 1
    K = full ^ H
    band = (n + 3) // 4
    allowed_left = 0
    allowed_right = 0
    for y in range(n):
        for x in range(n):
            v = y * n + x
            if (H >> v) & 1 and x < band and not (adj[v] & ~H):
                allowed_left |= 1 << v
            if (K >> v) & 1 and x >= n - band and not (adj[v] & ~K):
                allowed_right |= 1 << v
    return (h, H, K, allowed_left, allowed_right)


def event_info(word, data, n, adj):
    h, H, K, allowed_left, allowed_right = data
    if len(word) != n * n or len(set(word)) != n * n or set(word) != set(range(n * n)):
        return None
    pos = [0] * (n * n)
    for i, v in enumerate(word):
        pos[v] = i
    def unique_peak(G):
        peaks = []
        scan = G
        while scan:
            low = scan & -scan
            v = low.bit_length() - 1
            nbrs = adj[v] & G
            is_peak = True
            while nbrs:
                nbr_low = nbrs & -nbrs
                u = nbr_low.bit_length() - 1
                if pos[v] >= pos[u]:
                    is_peak = False
                    break
                nbrs ^= nbr_low
            if is_peak:
                peaks.append(v)
                if len(peaks) > 1:
                    return None
            scan ^= low
        return peaks[0] if len(peaks) == 1 else None
    lp = unique_peak(H)
    rp = unique_peak(K)
    if lp is None or rp is None or not ((allowed_left >> lp) & 1) or not ((allowed_right >> rp) & 1):
        return None
    # A global peak may be in each side only if both local peaks are roots.
    full = (1 << (n * n)) - 1
    gps = []
    for v in range(n * n):
        nbrs = adj[v]
        is_peak = True
        while nbrs:
            nbr_low = nbrs & -nbrs
            u = nbr_low.bit_length() - 1
            if pos[v] >= pos[u]:
                is_peak = False
                break
            nbrs ^= nbr_low
        if is_peak:
            gps.append(v)
    if len(gps) != 2:
        return None
    distance = abs(lp % n - rp % n) + abs(lp // n - rp // n)
    if distance < n - 2 * ((n + 3) // 4) + 1:
        return None
    return {"left_peak": lp, "right_peak": rp, "distance": distance, "curve": h}


def choose_chunks(curves, counts, total, master_seed, chunk_size):
    """Choose exactly total curve draws from integer C_h/S1, then group them."""
    S1 = sum(counts)
    rng = random.Random(master_seed)
    csum = []
    acc = 0
    for c in counts:
        acc += c
        csum.append(acc)
    selected = []
    for _ in range(total):
        r = rng.randrange(S1)
        lo, hi = 0, len(csum)
        while lo < hi:
            mid = (lo + hi) // 2
            if r < csum[mid]:
                hi = mid
            else:
                lo = mid + 1
        assert counts[lo] > 0
        selected.append(lo)
    # Curve-major chunks let a C++ counter cache serve at most 256 words.
    by_curve = defaultdict(int)
    for i in selected:
        by_curve[i] += 1
    chunks = []
    seq = 0
    for i, amount in sorted(by_curve.items()):
        while amount:
            k = min(chunk_size, amount)
            chunks.append({"job": seq, "curve_index": i, "curve": curves[i],
                           "samples": k, "cpp_seed": rng.getrandbits(64)})
            seq += 1
            amount -= k
    return chunks, Counter(selected)


def run_chunk(job, *, outdir, binary, n, all_curve_data, cache_mib, adj):
    raw = outdir / "raw" / f"job-{job['job']:05d}-{job['curve']}.json"
    stderr = outdir / "resource" / f"job-{job['job']:05d}-{job['curve']}.stderr"
    validation = outdir / "validation" / f"job-{job['job']:05d}-{job['curve']}.json"
    cmd = ["/usr/bin/time", "-l", str(binary), str(n), job["curve"], str(raw),
           "--samples", str(job["samples"]), "--seed", str(job["cpp_seed"]),
           "--cache-mib", str(cache_mib)]
    started = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    wrapper_pid = proc.pid
    stdout, stderr_text = proc.communicate()  # explicit wait/reap of each time+sampler child
    stderr.write_text(stderr_text)
    receipt = {**job, "command": cmd, "started_unix": started,
               "ended_unix": time.time(), "time_wrapper_pid": wrapper_pid,
               "waited_and_reaped": True, "returncode": proc.returncode,
               "resource_stderr_sha256": sha256(stderr)}
    if proc.returncode != 0 or not raw.exists():
        receipt["error"] = "sampler failed or did not emit raw JSON"
        return receipt, []
    doc = json.loads(raw.read_text())
    if int(doc.get("n", -1)) != n or doc.get("curve") != job["curve"] or len(doc.get("samples", [])) != job["samples"]:
        receipt["error"] = "sampler output schema/size mismatch"
        return receipt, []
    xs = []
    records = []
    for sample_index, word in enumerate(doc["samples"]):
        nominated = event_info(word, all_curve_data[job["curve_index"]], n, adj)
        if nominated is None:
            receipt["error"] = "nominated curve structural check failed"
            return receipt, []
        hits = []
        for i, d in enumerate(all_curve_data):
            if event_info(word, d, n, adj) is not None:
                hits.append(i)
        if job["curve_index"] not in hits or not hits:
            receipt["error"] = "chosen curve absent from exhaustive X count"
            return receipt, []
        xs.append(len(hits))
        # The word itself is retained in raw JSON.  This compact validation
        # record binds its index, all curve hits, and the nominated-event test.
        records.append({"sample_index": sample_index, "X": len(hits), "hits": hits,
                        "nominated": nominated})
    validation.write_text(json.dumps({"schema": "curved-mixture-exhaustive-X-v1", "job": job,
                                      "raw_sha256": sha256(raw), "samples": records},
                                     indent=2, sort_keys=True) + "\n")
    receipt.update({"raw_sha256": sha256(raw), "validation_sha256": sha256(validation),
                    "words_checked": len(xs), "status": "ok"})
    return receipt, xs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", type=Path)
    ap.add_argument("binary", type=Path)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--master-seed", type=int, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cache-mib", type=int, default=16)
    ap.add_argument("--chunk-size", type=int, default=256)
    args = ap.parse_args()
    if args.samples <= 0 or args.workers <= 0 or not 1 <= args.chunk_size <= 256 or args.cache_mib <= 0:
        raise SystemExit("invalid finite diagnostic parameters")
    if args.outdir.exists():
        raise SystemExit(f"refusing non-fresh output directory: {args.outdir}")
    if not args.inputs.is_file() or not args.binary.is_file():
        raise SystemExit("missing inputs or frozen binary")
    doc = json.loads(args.inputs.read_text())
    if doc.get("schema") != "exact-mixture-input-v1":
        raise SystemExit("unexpected exact input schema")
    n = int(doc["n"]); curves = list(doc["curves"]); counts = [int(x) for x in doc["counts"]]
    S1, S2, A2 = int(doc["S1"]), int(doc["S2"]), int(doc["A2"])
    if sum(counts) != S1 or len(curves) != len(counts) or any(c < 0 for c in counts):
        raise SystemExit("invalid curve count ledger")
    args.outdir.mkdir(parents=True)
    for d in ("raw", "resource", "validation", "receipts"):
        (args.outdir / d).mkdir()
    adj = adjacency(n)
    all_curve_data = [curve_data(n, h, adj) for h in curves]
    chunks, selected_counts = choose_chunks(curves, counts, args.samples, args.master_seed, args.chunk_size)
    manifest = {"schema": SCHEMA, "n": n, "inputs": str(args.inputs), "inputs_sha256": sha256(args.inputs),
                "binary": str(args.binary), "binary_sha256": sha256(args.binary),
                "samples": args.samples, "master_seed": args.master_seed, "workers": args.workers,
                "cache_mib": args.cache_mib, "chunk_size": args.chunk_size, "S1": str(S1), "S2": str(S2),
                "A2": str(A2), "selected_curve_counts": {curves[k]: v for k, v in selected_counts.items()},
                "jobs": chunks, "note": "finite PRNG diagnostic; no probability or asymptotic certification"}
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    all_info = []
    receipts = []
    # The exhaustive X pass is CPU-bound Python work; separate processes avoid
    # serializing it behind the GIL while C++ samplers run in their own children.
    with futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        pending = [pool.submit(run_chunk, j, outdir=args.outdir, binary=args.binary, n=n,
                               all_curve_data=all_curve_data, cache_mib=args.cache_mib, adj=adj) for j in chunks]
        for f in futures.as_completed(pending):
            rec, xs = f.result()
            receipts.append(rec)
            (args.outdir / "receipts" / f"job-{rec['job']:05d}.json").write_text(json.dumps(rec, indent=2, sort_keys=True)+"\n")
            if rec.get("status") != "ok":
                raise RuntimeError(json.dumps(rec, sort_keys=True))
            all_info.extend(xs)
    if len(all_info) != args.samples or sum(x <= 0 for x in all_info):
        raise RuntimeError("incomplete or zero-X aggregate")
    sum_x = sum(all_info)
    inv_sum = sum((Fraction(1, x) for x in all_info), Fraction(0))
    mean_inverse = inv_sum / args.samples
    estimator = Fraction(S1, A2) * mean_inverse
    hist = Counter(all_info)
    getcontext().prec = 50
    def dec(f): return str(Decimal(f.numerator) / Decimal(f.denominator))
    summary = {"schema": SCHEMA, "status": "complete", "n": n, "samples": args.samples,
               "S1": str(S1), "S2": str(S2), "A2": str(A2), "histogram_X": dict(sorted(hist.items())),
               "mean_X": {"fraction": f"{sum_x}/{args.samples}", "decimal": str(Decimal(sum_x)/Decimal(args.samples))},
               "known_S2_over_S1": (None if S2 == 0 else {"fraction": f"{S2}/{S1}", "decimal": dec(Fraction(S2,S1))}),
               "mean_inverse_X": {"fraction": f"{mean_inverse.numerator}/{mean_inverse.denominator}", "decimal": dec(mean_inverse)},
               "S1_over_A2_times_mean_inverse_X": {"fraction": f"{estimator.numerator}/{estimator.denominator}", "decimal": dec(estimator)},
               "receipt_count": len(receipts), "raw_count": len(list((args.outdir/'raw').glob('*.json'))),
               "validation_count": len(list((args.outdir/'validation').glob('*.json'))),
               "note": "Finite PRNG diagnostic only. The displayed ratio is not a rigorously certified probability estimate."}
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
