"""Verify distinct extended samples and summarize cross-machine steady measurements."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import statistics

from verify_and_aggregate import AFFINE_TOTAL, SUBSPACES, checked_shards, gauss, independent_truth, read_json


def batch_report(directory, label, weight_files, flat_files):
    weights = [read_json(path) for path in weight_files]
    flats = [read_json(path) for path in flat_files]
    intervals = [row for row in weights + flats]
    cpu_seconds = sum(row["cpu_seconds"] for row in intervals)
    busy_seconds = sum(row["wall_seconds"] for row in intervals)
    weight_cpu_seconds = sum(row["cpu_seconds"] for row in weights)
    weight_busy_seconds = sum(row["wall_seconds"] for row in weights)
    start = min(row["start_epoch"] for row in intervals)
    end = max(row["end_epoch"] for row in intervals)
    samples = []
    steady = []
    for file, row in zip(weight_files, weights):
        with file.with_name(file.stem + "-usage.csv").open() as stream:
            for data in csv.DictReader(stream):
                point = {k: float(v) for k, v in data.items()}
                samples.append(point)
                if 0.10 < point["elapsed_seconds"] < row["wall_seconds"] - 0.04:
                    steady.append(point)

    other_snapshots = []
    workloads = set()
    snapshot = []
    for line in (directory / ("extended_" + label) / "other_processes.txt").read_text().splitlines():
        fields = line.split(maxsplit=2)
        if len(fields) == 1:
            if snapshot:
                other_snapshots.append(snapshot)
            snapshot = []
            continue
        pname = Path(fields[2]).name
        if pname == "support_bitset":
            snapshot.append(float(fields[1]))
            workloads.add(fields[0])
    if snapshot:
        other_snapshots.append(snapshot)
    return {
        "workers": weights[0]["worker_count"], "logical_cpus": weights[0]["logical_cpu_count"],
        "start_epoch": start, "end_epoch": end, "elapsed_wall_seconds": end - start,
        "task_cpu_seconds": cpu_seconds, "task_wall_seconds": busy_seconds,
        "weight_task_cpu_seconds": weight_cpu_seconds, "weight_task_wall_seconds": weight_busy_seconds,
        "weight_average_task_core_equivalents": weight_cpu_seconds / weight_busy_seconds,
        "weight_p95_individual_cohort_core_equivalents": statistics.quantiles(
            [r["cpu_seconds"] / r["wall_seconds"] for r in weights], n=20
        )[18],
        "system_monitor_window_count": len(samples),
        "median_system_active_percent": statistics.median(p["system_active_percent"] for p in samples),
        "peak_task_core_equivalents_measured": max(p["task_core_equivalents"] for p in samples),
        "steady_window_count": len(steady),
        "steady_median_system_active_percent": (
            statistics.median(p["system_active_percent"] for p in steady) if steady else None
        ),
        "steady_median_task_core_equivalents": (
            statistics.median(p["task_core_equivalents"] for p in steady) if steady else None
        ),
        "other_support_workers_seen": len(workloads),
        "other_support_workers_peak_sum_reported_ps_percent": max(map(sum, other_snapshots), default=0),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    d = args.directory
    input_files = sorted((d / "extended_input").glob("cohort*.txt"))
    manifest = read_json(d / "extended_input" / "manifest.json")
    assert len(input_files) == manifest["batches"]
    labels = ("mini", "macbook", "mini_cross")
    weight_files = {h: [d / f"extended_{h}" / f"{p.stem}-weight.json" for p in input_files] for h in labels}
    flat_files = {h: [d / f"extended_{h}" / f"{p.stem}-flat.json" for p in input_files] for h in labels}
    rows = []
    duplicated_shards_checked = 0
    truth_inputs = set((d / "anf_samples.txt").read_text().splitlines())
    for j, p in enumerate(input_files):
        samples = p.read_text().splitlines()
        assert not truth_inputs.intersection(samples)
        truth_inputs.update(samples)
        z = {h: read_json(weight_files[h][j]) for h in labels}
        m = {h: read_json(flat_files[h][j]) for h in labels}
        wlist = checked_shards(z["mini"], z["macbook"])
        flist = checked_shards(m["mini"], m["macbook"])
        assert z["mini_cross"]["shards"] == z["macbook"]["shards"][:12]
        assert m["mini_cross"]["shards"] == m["macbook"]["shards"][:12]
        duplicated_shards_checked += 24
        mono = []
        tw = []
        for i, truth in enumerate(independent_truth(samples)):
            hist = [sum(part["histograms"][i][w] for part in wlist) for w in range(33)]
            assert sum(hist) == AFFINE_TOTAL and not any(hist[1::2])
            count = sum(part["samples"][i]["count"] for part in flist)
            assert hist[0] + hist[32] == count
            n = sum(truth)
            assert sum(w * m for w, m in enumerate(hist)) == n * SUBSPACES
            assert sum(w * (w - 1) * m for w, m in enumerate(hist)) == n * (n - 1) * gauss(8, 4)
            assert sum(w * (w - 1) * (w - 2) * m for w, m in enumerate(hist)) == n * (n - 1) * (n - 2) * gauss(7, 3)
            mono.append(count)
            tw.append(n)
        rows.append({
            "cohort": p.stem, "input_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "monochromatic_counts": mono, "truth_weights": tw,
        })
    assert len(truth_inputs) == manifest["samples"] + 64
    counts = Counter(a for row in rows for a in row["monochromatic_counts"])
    performance = {h: batch_report(d, h, weight_files[h], flat_files[h]) for h in labels}
    b, c = performance["macbook"], performance["mini_cross"]
    result = {
        **manifest, "affine_count_per_sample": AFFINE_TOTAL,
        "independent_histogram_and_bitplane_agreement": True,
        "checked_first_three_factorial_moments_for_all_samples": True,
        "duplicate_cross_machine_single_algorithm_shards_checked": duplicated_shards_checked,
        "monochromatic_count_distribution": dict(sorted(counts.items())),
        "performance": performance,
        "macbook_and_mini_cross_overlap_seconds": min(b["end_epoch"], c["end_epoch"]) - max(b["start_epoch"], c["start_epoch"]),
        "rows": rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
