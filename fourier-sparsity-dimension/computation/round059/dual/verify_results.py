"""Check archived receipts and independently audit selected covers over the rationals."""
import csv
import hashlib
import json
import statistics
import struct
import sys
from fractions import Fraction
from pathlib import Path

from make_bank import Generator, bank


KEYS = ("cases", "overlapping_cases", "t_equals_m_minus_maximum",
        "t_equals_centered_real_rank", "positional_checksum64", "by_dimension_and_family")
SOURCES = ("cover_cpu.c", "CoverMetal.swift", "bank.bin", "make_bank.py", "run_host.py")
MOD32 = (1 << 32) - 1
MOD64 = (1 << 64) - 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def hash32(value):
    value = ((value ^ (value >> 16)) * 0x7feb352d) & MOD32
    value = ((value ^ (value >> 15)) * 0x846ca68b) & MOD32
    return value ^ (value >> 16)


def hash64(case, details):
    value = (case << 24) | details
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MOD64
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MOD64
    return value ^ (value >> 31)


def generate(case, words):
    n = 3 + (case & 1 if case < 4 else (case >> 2) & 1)
    family = case & 3
    size = 1 << n
    full = (1 << size) - 1
    rng = Generator(hash32(case ^ 0xa2805901))
    weights = {}
    components = []

    def add(mask, weight):
        weights[mask] = weights.get(mask, 0) + weight

    def add_partition(index, weight):
        offset = 17 * (64 * (n - 3) + index)
        blocks = words[offset + 1:offset + 1 + words[offset]]
        assert sum(blocks) == full
        assert all(not (mask & other) for i, mask in enumerate(blocks) for other in blocks[i + 1:])
        for mask in blocks:
            add(mask, weight)

    if case < 2:
        target = sum(range(1, n + 1))
        for j in range(n):
            mask = sum(1 << x for x in range(size) if x >> j & 1)
            add(mask, j + 1)
            add(full ^ mask, j + 1)
    elif case < 4 or family >= 2:
        points = list(range(size))
        for i in range(size - 1, 0, -1):
            j = rng.get() % (i + 1)
            points[i], points[j] = points[j], points[i]
        target = 2
        start = 0
        for length in [3, 5] * (size // 8):
            cycle = points[start:start + length]
            components.append(cycle)
            for a, b in zip(cycle, cycle[1:] + cycle[:1]):
                add((1 << a) | (1 << b), 1)
            start += length
        if case >= 4 and family == 3:
            add_partition(rng.get() % 64, 2)
            target += 2
    else:
        target = 0
        for _ in range(2 + family):
            index, weight = rng.get() % 64, rng.get() % 5 + 1
            add_partition(index, weight)
            target += weight
    return n, weights, target, components


def q_rank(matrix):
    basis = {}
    for row in matrix:
        residue = list(map(Fraction, row))
        for index, b in sorted(basis.items()):
            coefficient = residue[index]
            if coefficient:
                residue = [v - coefficient * w for v, w in zip(residue, b)]
        if any(residue):
            index = next(i for i, v in enumerate(residue) if v)
            basis[index] = [v / residue[index] for v in residue]
    return len(basis)


def audit_case(case, words):
    n, weights, target, components = generate(case, words)
    N = 1 << n
    masks = list(weights)
    incidence = [[int(mask >> x & 1) for x in range(N)] for mask in masks]
    assert all(weight > 0 for weight in weights.values())
    assert all(sum(value * weights[mask] for value, mask in zip(col, masks)) == target
               for col in zip(*incidence))
    normal_span = {0}
    for row in incidence:
        points = [i for i, value in enumerate(row) if value]
        assert all(x ^ y ^ points[0] in points for x in points for y in points)
        for a in range(1, N):
            if len({(a & point).bit_count() % 2 for point in points}) == 1:
                normal_span |= {a ^ g for g in normal_span}
    t = len(normal_span).bit_length() - 1
    rank = q_rank([[a - row[0] for a in row[1:]] for row in incidence])
    max_active = max(map(sum, zip(*incidence)))
    m = len(masks)
    assert t <= min(m - max_active, rank)
    code = m | max_active << 6 | t << 12 | rank << 15 | (n - 3) << 20 | (case & 3) << 21
    cert = dict(id=case, n=n, m=m, max_active=max_active, t=t, real_centered_rank=rank,
                mask_weight_pairs=list(weights.items()), weights_divided_by=target, packed_code=code)
    if case in (2, 3):
        attainable = {0}
        for mask in masks:
            assert mask.bit_count() == 2
            attainable |= {mask | prev for prev in attainable if not (mask & prev)}
        assert (1 << N) - 1 not in attainable
        assert all(len(cycle) % 2 == 1 for cycle in components)
        cert.update(odd_components=components, partition_subset_exists=False)
    return cert


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        assert a <= b
        if result and a <= result[-1][1]:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return result


def intersect(left, right):
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            result.append([a, b])
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return result


def summarize(path, source_hashes, epoch):
    cpu, gpu, receipt = [read_json(path / x) for x in ("cpu.json", "gpu.json", "receipt.json")]
    assert receipt["cpu_gpu_exit_codes"] == [0, 0] and not receipt["remaining_task_process_pids"]
    assert receipt["source_sha256"] == source_hashes
    assert receipt["same_requested_epoch"] == epoch
    assert receipt["cpu_gpu_outputs_sha256"] == [sha(path / x) for x in ("cpu.json", "gpu.json")]
    assert receipt["sampled_cpu_telemetry_sha256"] == sha(path / "cpu_usage.csv")
    first, last = receipt["input_interval_half_open"]
    assert cpu["cases"] == gpu["cases"] == cpu["count"] == gpu["count"] == last - first
    assert cpu["first_id"] == gpu["first_id"] == first
    assert all(cpu[key] == gpu[key] for key in KEYS)
    assert sum(map(sum, cpu["by_dimension_and_family"])) == last - first
    rows = list(csv.DictReader((path / "cpu_usage.csv").read_text().splitlines()))
    full_cpu_samples = [x for x in rows if int(x["active_workers"]) == cpu["threads"]]
    assert cpu["threads"] == receipt["logical_cores"] == len(cpu["worker_cpu_seconds"])
    worker_utilization = [float(x["a28_process_cpu_percent"]) for x in full_cpu_samples]
    machine_utilization = [float(x["whole_machine_cpu_active_percent"]) for x in full_cpu_samples]
    gpu_utilization = [n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                       for n in sample["metrics"]["Device Utilization %"]]
    gpu_intervals = union(gpu["gpu_command_intervals_epoch"])
    assert len(gpu["gpu_command_intervals_epoch"]) == gpu["gpu_command_buffers"]
    assert gpu["all_checks_passed"]
    summary = {key: cpu[key] for key in KEYS}
    summary.update(interval_half_open=[first,last], logical_cores=cpu["threads"],
                   a28_cpu_avg_core_equivalents=cpu["cpu_seconds"]/cpu["wall"],
                   all_worker_cpu_seconds=[min(cpu["worker_cpu_seconds"]),max(cpu["worker_cpu_seconds"])],
                   a28_cpu_measured_peak_percent=max(worker_utilization, default=None),
                   a28_cpu_measured_median_percent=statistics.median(worker_utilization) if worker_utilization else None,
                   machine_cpu_measured_peak_percent=max(machine_utilization, default=None),
                   full_worker_cpu_sample_count=len(full_cpu_samples),
                   device=gpu["device"], gpu_command_buffers=gpu["gpu_command_buffers"],
                   gpu_command_seconds=gpu["gpu_command_seconds"], gpu_wall_seconds=gpu["metal_wall_seconds"],
                   gpu_device_sample_count=len(gpu_utilization),
                   gpu_device_peak_percent=max(gpu_utilization, default=None),
                   gpu_device_median_percent=statistics.median(gpu_utilization) if gpu_utilization else None,
                   gpu_baseline_percent=receipt["GPU_before"]["metrics"]["Device Utilization %"],
                   cpu_pid=receipt["cpu_pid"], metal_pid=receipt["metal_pid"],
                   exit_codes=receipt["cpu_gpu_exit_codes"], remaining_pids=receipt["remaining_task_process_pids"])
    return summary, [[cpu["start_epoch"], cpu["end_epoch"]]], gpu_intervals, gpu["witness_codes"]


def main(local, macbook, output):
    root = Path(__file__).resolve().parent
    words = bank()
    assert (root / "bank.bin").read_bytes() == struct.pack("<%dI" % len(words), *words)
    hashes = {name: sha(root / name) for name in SOURCES}
    start = read_json(root / "START.json")
    epoch = start["common_epoch"]
    mini, cpu_a, gpu_a, witness_a = summarize(local, hashes, epoch)
    laptop, cpu_b, gpu_b, witness_b = summarize(macbook, hashes, epoch)
    assert mini["interval_half_open"][0] == 0
    assert mini["interval_half_open"][1] == laptop["interval_half_open"][0]
    assert not witness_b
    designated = [audit_case(i, words) for i in range(4)]
    assert {str(item["id"]): item["packed_code"] for item in designated} == witness_a
    sampled = [audit_case(i, words) for i in range(16384)]
    pilot_checksum = 0
    for sample in sampled:
        pilot_checksum ^= hash64(sample["id"], sample["packed_code"])
    pilot = read_json(local / "pilot_cpu.json")
    assert all(pilot[key] == read_json(local / "pilot_gpu.json")[key] for key in KEYS)
    assert f"{pilot_checksum:x}" == pilot["positional_checksum64"]
    assert sum(x["max_active"] > 1 for x in sampled) == pilot["overlapping_cases"]
    simultaneous = intersect(intersect(cpu_a, cpu_b), intersect(gpu_a, gpu_b))
    output.write_text(json.dumps(dict(common_epoch=epoch, local=mini, macbook=laptop,
                                     four_way_active_seconds=sum(b - a for a, b in simultaneous),
                                     four_way_positive_intervals=len(simultaneous),
                                     independent_fraction_cases=len(sampled),
                                     designated_covers=designated), indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main(*(Path(item) for item in sys.argv[1:]))
