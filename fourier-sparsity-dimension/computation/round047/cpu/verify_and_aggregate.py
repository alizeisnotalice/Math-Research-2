"""Check both enumerations, witnesses, and exact incidence moments."""
import argparse
from collections import Counter
import csv
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import prod
from pathlib import Path
import statistics


N = 9
D = 5
X = 1 << N
FIBRE = 1 << (N - D)
TERMS = [i for i in range(X) if i.bit_count() <= 3]


@lru_cache
def gauss(n, k):
    if not 0 <= k <= n:
        return 0
    if not k:
        return 1
    return gauss(n - 1, k) + (1 << (n - k)) * gauss(n - 1, k - 1)


SUBSPACES = gauss(N, D)
AFFINE_TOTAL = SUBSPACES * FIBRE
assert SUBSPACES == prod(X - (1 << i) for i in range(D)) // prod((1 << D) - (1 << i) for i in range(D))


def read_json(path):
    return json.loads(path.read_bytes())


def independent_truth(samples):
    assert len(samples) == 64 and len(set(samples)) == 64
    result = []
    for row in samples:
        assert len(row) == len(TERMS) and set(row) <= {"0", "1"}
        terms = [term for term, coefficient in zip(TERMS, row) if coefficient == "1"]
        assert max(term.bit_count() for term in terms) == 3
        result.append([sum(x & term == term for term in terms) % 2 for x in range(X)])
    return result


def pattern_index():
    index = []
    pos = 0
    for pivot_mask in sorted(sum(1 << i for i in pivot) for pivot in combinations(range(N), D)):
        pivot = [i for i in range(N) if (pivot_mask >> i) & 1]
        free = [i for i in range(N) if (pivot_mask >> i) & 1 == 0]
        options = [[j for j in free if j > i] for i in pivot]
        new_pos = pos + (1 << sum(map(len, options)))
        index.append((pos, new_pos, pivot, free, options))
        pos = new_pos
    assert pos == SUBSPACES
    return index


INDEX = pattern_index()


def decode(flat_id):
    ordinal, fibre = divmod(flat_id, FIBRE)
    lo, _, pivot, free, options = next(item for item in INDEX if item[0] <= ordinal < item[1])
    digit = ordinal - lo
    rows = []
    for p, columns in zip(pivot, options):
        digit, local = divmod(digit, 1 << len(columns))
        row = (1 << p) | sum(1 << a for j, a in enumerate(columns) if (local >> j) & 1 == 0)
        rows.append(row)
    assert not digit
    offset = sum(1 << a for j, a in enumerate(free) if (fibre >> j) & 1 == 0)
    return rows, offset


def points(rows, offset):
    vertices = {offset}
    for r in rows:
        vertices |= {v ^ r for v in vertices}
    assert len(vertices) == 32
    return vertices


def checked_shards(a, b):
    parts = sorted(a["shards"] + b["shards"], key=lambda shard: shard["id"])
    assert [p["id"] for p in parts] == list(range(27))
    for i, part in enumerate(parts):
        lo, hi = SUBSPACES * i // 27, SUBSPACES * (i + 1) // 27
        assert (part["linear_from"], part["linear_to"], part["affine_count"]) == (lo, hi, (hi - lo) * FIBRE)
    assert sum(p["affine_count"] for p in parts) == AFFINE_TOTAL
    return parts


def monitor_summary(report, usage):
    with usage.open() as stream:
        samples = list(csv.DictReader(stream))
    active = [float(row["system_active_percent"]) for row in samples]
    task = [float(row["task_core_equivalents"]) for row in samples]
    all_elapsed = [float(row["elapsed_seconds"]) for row in samples]
    middle = [i for i, elapsed in enumerate(all_elapsed) if 0.10 < elapsed < report["wall_seconds"] - 0.04]
    return {
        "workers": report["worker_count"], "cores": report["logical_cpu_count"],
        "seconds": report["wall_seconds"], "cpu_seconds": report["cpu_seconds"],
        "average_task_core_equivalents": report["cpu_seconds"] / report["wall_seconds"],
        "monitor_samples": len(samples), "system_active_percent_median": statistics.median(active),
        "task_core_equivalents_peak_window": max(task), "steady_window_samples": len(middle),
        "steady_system_active_percent_median": statistics.median(active[i] for i in middle) if middle else None,
        "steady_task_core_equivalents_median": statistics.median(task[i] for i in middle) if middle else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    loc = args.directory
    mini, mac = [read_json(loc / f"{host}_shards.json") for host in ("mini", "macbook")]
    mini_weights, mac_weights = [read_json(loc / f"{host}_weights.json") for host in ("mini", "macbook")]
    fast = checked_shards(mini, mac)
    full = checked_shards(mini_weights, mac_weights)
    s_path = loc / "anf_samples.txt"
    truth = independent_truth(s_path.read_text().splitlines())
    pair_incidence = gauss(8, 4)
    plane_incidence = gauss(7, 3)
    triple_dim_incidence = gauss(6, 2)
    pilots = read_json(loc / "pilot.json")["shards"][0]
    bitplane = [sum(t[x] << j for j, t in enumerate(truth)) for x in range(X)]
    expected_pilot = [0] * 64
    for flat_id in range(FIBRE * pilots["linear_from"], FIBRE * pilots["linear_to"]):
        rows, offset = decode(flat_id)
        any_one, all_one = 0, (1 << 64) - 1
        for x in points(rows, offset):
            any_one |= bitplane[x]
            all_one &= bitplane[x]
        hits = all_one | ((1 << 64) - 1) ^ any_one
        while hits:
            bit = hits & -hits
            expected_pilot[bit.bit_length() - 1] += 1
            hits ^= bit
    assert expected_pilot == [r["count"] for r in pilots["samples"]]

    rows = []
    verified_witnesses = 0
    for i, bits in enumerate(truth):
        hist = [sum(shard["histograms"][i][j] for shard in full) for j in range(33)]
        assert sum(hist) == AFFINE_TOTAL and all(hist[j] == 0 for j in range(1, 33, 2))
        support = [x for x, bit in enumerate(bits) if bit]
        support_size = len(support)
        f1 = sum(k * n for k, n in enumerate(hist))
        f2 = sum(k * (k - 1) * n for k, n in enumerate(hist))
        f3 = sum(k * (k - 1) * (k - 2) * n for k, n in enumerate(hist))
        f4 = sum(k * (k - 1) * (k - 2) * (k - 3) * n for k, n in enumerate(hist))
        assert f1 == support_size * SUBSPACES
        assert f2 == support_size * (support_size - 1) * pair_incidence
        assert f3 == support_size * (support_size - 1) * (support_size - 2) * plane_incidence
        dependencies = 0
        for delta in range(1, X):
            intersection = sum(bits[x ^ delta] for x in support)
            dependencies += intersection * (intersection - 2)
        ordered_distinct = prod(support_size - j for j in range(4))
        assert f4 == dependencies * plane_incidence + (ordered_distinct - dependencies) * triple_dim_incidence

        count = 0
        witnesses = []
        for shard in fast:
            entry = shard["samples"][i]
            count += entry["count"]
            w = entry["witness"]
            if not entry["count"]:
                assert w is None
                continue
            witness_id = w["flat_id"]
            assert FIBRE * shard["linear_from"] <= witness_id < FIBRE * shard["linear_to"]
            matrix, offset = decode(witness_id)
            assert w["rows"] == matrix and w["offset"] == offset
            assert {bits[x] for x in points(matrix, offset)} == {w["bit"]}
            verified_witnesses += 1
            witnesses.append(w)
        assert count == hist[0] + hist[-1]
        rows.append({
            "sample": i, "truth_weight": support_size, "monochromatic_flat_count": count,
            "zero_flats": hist[0], "one_flats": hist[-1], "restriction_weight_histogram": hist,
            "checked_incidence_factorial_moments_1_to_4": [f1, f2, f3, f4],
            "first_witness": min(witnesses, key=lambda w: w["flat_id"], default=None),
        })
    results = {
        "input_sha256": hashlib.sha256(s_path.read_bytes()).hexdigest(),
        "linear_subspaces_gaussian": SUBSPACES, "affine_flats": AFFINE_TOTAL,
        "checked_pilot_affine_flats": pilots["affine_count"],
        "verified_positive_shard_witnesses": verified_witnesses,
        "sample_mono_count_distribution": dict(sorted(Counter(row["monochromatic_flat_count"] for row in rows).items())),
        "sample_ids_with_no_monochromatic_5flat": [row["sample"] for row in rows if row["monochromatic_flat_count"] == 0],
        "runs": {
            "first_mini": monitor_summary(mini, loc / "mini_usage.csv"),
            "first_macbook": monitor_summary(mac, loc / "macbook_usage.csv"),
            "weight_mini": monitor_summary(mini_weights, loc / "mini_weights_usage.csv"),
            "weight_macbook": monitor_summary(mac_weights, loc / "macbook_weights_usage.csv"),
        },
        "weight_run_overlap_seconds": min(mini_weights["end_epoch"], mac_weights["end_epoch"])
        - max(mini_weights["start_epoch"], mac_weights["start_epoch"]),
        "rows": rows,
    }
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
