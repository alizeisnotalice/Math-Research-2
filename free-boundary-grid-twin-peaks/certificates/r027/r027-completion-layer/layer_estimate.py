#!/usr/bin/env python3
"""Exact row-transfer layer census and conservative memory arithmetic.

No completion recurrence is evaluated here.  The forbidden/IE set after the
top-row prefix in P_n is an (n-2)-by-n rectangle.
"""
from __future__ import annotations

from collections import defaultdict
from math import factorial, ceil, log2
from pathlib import Path
import json


def layer_hist(rows: int, cols: int) -> dict[int, int]:
    valid = [s for s in range(1 << cols) if not (s & (s << 1))]
    state = {(0, 0): 1}
    for _ in range(rows):
        nxt = defaultdict(int)
        for (previous, k), count in state.items():
            for row in valid:
                if not (previous & row):
                    nxt[row, k + row.bit_count()] += count
        state = nxt
    out = defaultdict(int)
    for (_, k), count in state.items():
        out[k] += count
    return dict(out)


def fmt(n: int) -> dict:
    return {"bytes": n, "MiB": n / 2**20, "GiB": n / 2**30}


def report(n: int) -> dict:
    h = layer_hist(n - 2, n)
    # A C(I) is U256=32 B.  Rank-addressed arrays store values only; a
    # compact open-addressed predecessor index needs a 64-bit key + 32-bit
    # index but rounds to 16 bytes/slot.  Capacity is ceil(L/0.70).
    pairs = []
    for k in range(1, max(h) + 1):
        prev, cur = h.get(k - 1, 0), h[k]
        cap = ceil(prev / 0.70)
        pairs.append({"k": k, "previous": prev, "current": cur,
                      "two_value_arrays": 32 * (prev + cur),
                      "compact_hash_previous": 16 * cap,
                      "hash_plus_two_arrays": 32 * (prev + cur) + 16 * cap,
                      "sorted_masks_plus_two_values": 8 * prev + 32 * (prev + cur)})
    peak_hash = max(pairs, key=lambda x: x["hash_plus_two_arrays"])
    peak_sorted = max(pairs, key=lambda x: x["sorted_masks_plus_two_values"])
    predecessor_lookups = sum(k * h[k] for k in h)
    binary_comparisons_upper = sum(k * h[k] * max(1, (h.get(k - 1, 1) - 1).bit_length()) for k in range(1, max(h) + 1))
    # A power-of-two, parallel-array open-addressed table with 64-bit key and
    # 32-bit value uses 12 B/slot and avoids C++ struct padding.  Rebuilding
    # the table by a deterministic second enumeration of the layer avoids a
    # current-layer mask array.
    p = peak_hash["previous"]
    pow2cap = 1 << (p - 1).bit_length()
    parallel_hash_bytes = 12 * pow2cap
    reenumerated_hash_peak = peak_hash["two_value_arrays"] + parallel_hash_bytes
    return {"n": n, "rectangle": [n - 2, n],
            "valid_row_masks": len([s for s in range(1 << n) if not (s & (s << 1))]),
            "layer_histogram": {str(k): h[k] for k in sorted(h)},
            "total_independent_sets": sum(h.values()),
            "max_layer": max(h, key=h.get), "max_layer_count": max(h.values()),
            "u256_C_bits": 256, "C_upper_bits": factorial(n * (n - 1)).bit_length(),
            "coarse_signed_sum_bits": (sum(h.values()) * factorial(n * (n - 1))).bit_length(),
            "peak_compact_hash_model": {**peak_hash, "formatted": fmt(peak_hash["hash_plus_two_arrays"])},
            "peak_sorted_binary_model": {**peak_sorted, "formatted": fmt(peak_sorted["sorted_masks_plus_two_values"])},
            "predecessor_lookups": predecessor_lookups,
            "sorted_binary_comparisons_upper": binary_comparisons_upper,
            "re_enumerated_parallel_hash_model": {"transition_k": peak_hash["k"], "previous": p,
                "power_of_two_slots": pow2cap, "load_factor": p / pow2cap,
                "two_value_arrays": peak_hash["two_value_arrays"], "parallel_hash_bytes": parallel_hash_bytes,
                "total": reenumerated_hash_peak, "formatted": fmt(reenumerated_hash_peak)},
            "layer_pairs": pairs}


out = {"method": "row-transfer exact layer census; memory models exclude allocator/OS/thread overhead",
       "n3_to_n8": [report(n) for n in range(3, 9)]}
Path("layer-estimate.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
