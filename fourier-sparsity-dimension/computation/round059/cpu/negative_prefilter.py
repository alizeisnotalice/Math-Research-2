"""Exhaust the Fourier sign obstruction for a cover missing few domains."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import time
from check import affine_flats, binary_rank


def profile(n):
    flats = affine_flats(n)
    active_flats = []
    missed_flats = []
    for points, normals in flats:
        mask = sum(1 << w for w in normals)
        if 0 in points:
            active_flats.append((points, mask))
        else:
            negative = sum(1 << w for w in normals if (w & points[0]).bit_count() % 2)
            missed_flats.append((points, mask, negative))
    ranks = {}
    counts = Counter()
    first_possible = {}
    for k in range(1, n):
        for flats in combinations(missed_flats, k):
            normals = 0
            negatives = 0
            for _, nm, neg in flats:
                normals |= nm
                negatives |= neg
            if normals != negatives:
                continue
            if normals not in ranks:
                ranks[normals] = binary_rank(w for w in range(1, 1 << n) if (normals >> w) & 1)
            rank = ranks[normals]
            counts[(k, rank)] += 1
            if rank > k and k not in first_possible:
                first_possible[k] = {"rank": rank, "missed_flats": [points for points, _, _ in flats]}
    return {"n": n, "active_flats": len(active_flats), "missed_flats": len(missed_flats),
            "compatible_negative_candidate_counts": {f"{k},{t}": c for (k,t),c in sorted(counts.items())},
            "first_possible_violations": first_possible}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    result = [profile(n) for n in range(1, 5)]
    payload = {"census": result, "elapsed_seconds": time.perf_counter() - start}
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(payload)
