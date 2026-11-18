"""Independently reproduce Walsh checks and exact restriction expectations."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path


KEYS = ("first_id", "count", "rank_losses", "cases_with_collisions", "cases_with_missing_basis_bands",
        "sums", "positional_checksum64", "tested_by_hk", "collision_pairs_by_hk")
SOURCES = ("restriction_cpu.c", "RestrictionCommon.h", "RestrictionKernel.metal", "RestrictionMetal.swift",
           "run_host.py")
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def mix(x):
    x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & MASK64
    return x ^ (x >> 31)


def parity(x):
    return x.bit_count() & 1


def characters(mask, arg):
    return 1 - 2 * parity(mask & arg)


def span(vectors):
    result = {0}
    for v in vectors:
        result.update([u ^ v for u in result])
    return result


def model(case):
    h = 5 + (case >> 1 & 1)
    k = 1 + (case >> 2) % (h - 1)
    holes = case >> 3 & 1
    if case < 4:
        return h, [2, 4], case, [2, 4 | 1 << h, 8 | 2 << h], 0x97, holes
    state = case ^ 0xa2806301 or 1

    def word():
        nonlocal state
        state ^= state << 13
        state ^= (state & MASK32) >> 17
        state ^= state << 5
        state &= MASK32
        return state

    mask = (1 << h) - 1
    shift = word() & mask
    forms = [word() & (mask ^ 1) | (i << h) for i in range(3)]
    table = word() & 255
    if holes:
        table &= 127
        table |= parity(table) << 7
    basis = []
    while len(basis) < k:
        v = word() & mask
        if v not in span(basis):
            basis.append(v)
    return h, basis, shift, forms, table, holes


def truth(h, forms, table, holes):
    points = []
    for x in range(8 << h):
        if x & 1:
            idx = sum(parity(form & x) << j for j, form in enumerate(forms))
            negative = ((table >> idx) ^ (x >> (h + 2))) & 1
        else:
            y = x & ((1 << h) - 2)
            negative = y in (0, 2, 4, 8) if holes else y == 0
        points.append(1 - 2 * negative)
    return points


def transform(values):
    if len(values) == 1:
        return values
    first = transform(values[::2])
    second = transform(values[1::2])
    answer = []
    for a, b in zip(first, second):
        answer.extend([a + b, a - b])
    return answer


def restricted_truth(values, h, rows, shift):
    sites = [shift]
    for v in rows:
        sites.extend([s ^ v for s in sites])
    return [values[(beta << h) | s] for beta in range(8) for s in sites]


def view_spectrum(spectrum, n):
    result = [[] for _ in range(8)]
    for beta in range(8):
        result[beta] = [a for a in range(1 << n) if spectrum[(beta << n) | a]]
    return result


def representative_basis(bands):
    choices = []
    for beta, band in enumerate(bands):
        if beta and band and beta not in span([item[0] for item in choices]):
            choices.append((beta, band[0]))
    return choices


def maps(h, rows):
    return [sum(parity(a & r) << j for j, r in enumerate(rows)) for a in range(1 << h)]


def audit(case):
    h, rows, shift, forms, table, holes = model(case)
    k = len(rows)
    M, V = 1 << h, 1 << k
    original_truth = truth(h, forms, table, holes)
    w = transform(original_truth)
    actual = transform(restricted_truth(original_truth, h, rows, shift))
    image = maps(h, rows)
    folded = [sum(w[(beta << h) | a] * characters(a, shift) for a in range(M) if image[a] == output)
              for beta in range(8) for output in range(V)]
    assert [v * (M // V) for v in actual] == folded
    step = 16 << holes
    assert all(c % step == 0 for c in w)
    core = [i for i, c in enumerate(w[:M]) if c]
    odd = [i for i, c in enumerate(w) if (c // step) % 2]
    if holes:
        assert core == [i for i in range(M) if (i & 14).bit_count() != 2]
        assert odd == [i for i in range(M) if parity(i & 14)]
    else:
        assert core == odd == list(range(M))
    assert len(span(core)) == M
    old = view_spectrum(w, h)
    new = view_spectrum(actual, k)
    basis = representative_basis(old)
    eligible = sum(len(old[beta]) - 1 for beta, _ in basis)
    pair_counts = [sum(a != rep and image[a] == image[rep] for a in old[beta]) for beta, rep in basis]
    missing = sum(not new[beta] for beta, _ in basis)
    matched = sum(c != 0 for c in pair_counts)
    original = sum(map(len, old))
    before, after = len(basis), len(representative_basis(new))
    assert 0 <= before - after <= missing <= matched <= sum(pair_counts) <= eligible <= original - len(core) - before
    spectral = 0
    for v in actual:
        spectral = (((spectral * 16777619) & MASK32) ^ (v & MASK32))
    statistics_row = [original, sum(map(len, new)), before, after, eligible, sum(pair_counts),
                      matched, missing, len(core), h, k, spectral, len(odd), holes]
    certificate = dict(id=case, h=h, k=k, n=h + 3, shift=shift, rows=rows,
                       phase_forms=forms, phase_table=table, holes=bool(holes), spectral_core=core, highest_odd_layer=odd,
                       external_basis_representatives=basis,
                       expected_diagnostics=statistics_row)
    return statistics_row, certificate


def subspaces(h, k):
    for leading in combinations(range(h), k):
        others = [j for j in range(h) if j not in leading]
        parameters = [(row, j) for row in range(k) for j in others if j > leading[row]]
        for free_bits in product((0, 1), repeat=len(parameters)):
            rows = [1 << i for i in leading]
            for (row, j), value in zip(parameters, free_bits):
                rows[row] |= value << j
            yield rows


def coset_representatives(h, rows):
    count = 1 << h
    unused = set(range(count))
    space = span(rows)
    while unused:
        a = min(unused)
        yield a
        unused.difference_update(a ^ point for point in space)


def expected(case):
    h, _, _, forms, table, holes = model(case)
    tt = truth(h, forms, table, holes)
    old = view_spectrum(transform(tt), h)
    basis = representative_basis(old)
    d = len(basis)
    eligible = sum(len(old[beta]) - 1 for beta, _ in basis)
    values = []
    for k in range(1, h):
        geometry_count = [0] * (1 << h)
        distribution = Counter()
        number_w = number_affine = same_representative_pairs = 0
        for rows in subspaces(h, k):
            number_w += 1
            projected = maps(h, rows)
            for diff in range(1, 1 << h):
                geometry_count[diff] += (projected[diff] == 0)
            pairs = sum(a != rep and projected[a] == projected[rep]
                        for beta, rep in basis for a in old[beta])
            colliding = sum(any(a != rep and projected[a] == projected[rep] for a in old[beta])
                            for beta, rep in basis)
            same_representative_pairs += pairs
            for shift in coset_representatives(h, rows):
                number_affine += 1
                actual = view_spectrum(transform(restricted_truth(tt, h, rows, shift)), k)
                lost = d - len(representative_basis(actual))
                assert 0 <= lost <= colliding <= pairs
                distribution[lost] += 1
        ratio = Fraction((1 << (h-k))-1, (1 << h)-1)
        assert len(set(geometry_count[1:])) == 1 and Fraction(geometry_count[1],number_w) == ratio
        assert number_affine == number_w * (1 << (h-k))
        mean_pairs = Fraction(same_representative_pairs, number_w)
        assert mean_pairs == eligible * ratio
        mean_loss = Fraction(sum(key * v for key, v in distribution.items()),number_affine)
        assert mean_loss <= mean_pairs
        values.append(dict(k=k, all_linear_subspaces=number_w, all_cosets=number_affine,
                           nonzero_difference_collision_count=geometry_count[1],
                           one_pair_collision_probability=str(ratio), eligible_pairs=eligible,
                           expected_representative_collision_pairs=str(mean_pairs),
                           actual_rank_loss_histogram=dict(sorted(distribution.items())),
                           expected_actual_external_rank_loss=str(mean_loss)))
    return dict(h=h, model_id=case, core_support_size=len(old[0]), basis=basis, table=table, forms=forms, rows=values)


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        assert a <= b
        if result and a <= result[-1][1]:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return result


def intersect(a, b):
    result = []
    i = j = 0
    while i < len(a) and j < len(b):
        start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if start < end:
            result.append([start, end])
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return result


def receipt_summary(path, hashes, epoch):
    cpu, gpu, receipt = [read(path / x) for x in ("cpu.json", "gpu.json", "receipt.json")]
    assert all(cpu[key] == gpu[key] for key in KEYS)
    assert receipt["same_requested_epoch"] == epoch
    assert receipt["source_sha256"] == hashes
    assert receipt["cpu_gpu_exit_codes"] == [0,0] and not receipt["remaining_task_process_pids"]
    assert receipt["cpu_gpu_outputs_sha256"] == [sha(path / x) for x in ("cpu.json", "gpu.json")]
    assert receipt["sampled_cpu_telemetry_sha256"] == sha(path / "cpu_usage.csv")
    first, last = receipt["input_interval_half_open"]
    assert first == cpu["first_id"] and last-first == cpu["count"]
    assert sum(map(sum, cpu["tested_by_hk"])) == cpu["count"]
    assert sum(map(sum, cpu["collision_pairs_by_hk"])) == cpu["sums"][5]
    assert all(v > 0 for v in cpu["worker_cpu_seconds"])
    assert cpu["threads"] == receipt["logical_cores"] == len(cpu["worker_cpu_seconds"])
    assert gpu["all_checks_passed"] and len(gpu["gpu_command_intervals_epoch"]) == gpu["gpu_command_buffers"]
    csv_rows = list(csv.DictReader((path / "cpu_usage.csv").read_text().splitlines()))
    active = [row for row in csv_rows if int(row["active_workers"]) == cpu["threads"]]
    cpu_load = [float(r["a28_process_cpu_percent"]) for r in active]
    idle = [100.0-float(r["whole_machine_cpu_active_percent"]) for r in active
            if float(r["whole_machine_cpu_active_percent"]) >= 0]
    gpu_load = [n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                for n in sample["metrics"]["Device Utilization %"]]
    names = ("s_original", "s_restricted", "external_rank_before", "external_rank_after",
             "eligible_basis_pairs", "colliding_basis_pairs", "collided_basis_representatives",
             "missing_basis_bands", "core_support_size", "h", "k", "compressed_spectral_fingerprint",
             "highest_odd_layer_size", "holed_case_count")
    summary = {key: cpu[key] for key in KEYS if key != "sums"}
    summary.update(sum_diagnostics=dict(zip(names,cpu["sums"])),
                   actual_external_rank_loss_sum=cpu["sums"][2]-cpu["sums"][3],
                   interval_half_open=[first,last], cores=cpu["threads"],
                   average_a28_cpu_core_equivalents=cpu["cpu_seconds"]/cpu["wall"],
                   peak_a28_cpu_percent=max(cpu_load,default=None),
                   median_a28_cpu_percent=statistics.median(cpu_load) if cpu_load else None,
                   full_worker_cpu_samples=len(cpu_load),
                   worker_wall_seconds=cpu["wall"], total_project_cpu_seconds=cpu["cpu_seconds"],
                   whole_machine_idle_sample_count=len(idle),
                   whole_machine_idle_median_percent=statistics.median(idle) if idle else None,
                   whole_machine_idle_min_percent=min(idle, default=None),
                   whole_machine_idle_max_percent=max(idle, default=None),
                   zero_idle_samples=sum(v==0 for v in idle),
                   peak_whole_machine_cpu_percent=receipt.get("peak_whole_machine_cpu_percent"),
                   device=gpu["device"], gpu_command_count=gpu["gpu_command_buffers"],
                   gpu_command_duration_sum=gpu["gpu_command_seconds"],
                   metal_wall_seconds=gpu["metal_wall_seconds"], gpu_device_samples=len(gpu_load),
                   gpu_device_peak_percent=max(gpu_load,default=None),
                   gpu_device_median_percent=statistics.median(gpu_load) if gpu_load else None,
                   gpu_device_baseline=receipt["GPU_before"]["metrics"]["Device Utilization %"],
                   exit_codes=receipt["cpu_gpu_exit_codes"], remaining_pids=receipt["remaining_task_process_pids"])
    return summary, [[cpu["start_epoch"],cpu["end_epoch"]]], union(gpu["gpu_command_intervals_epoch"]),gpu["controls"]


def main(local, macbook, output):
    root=Path(__file__).resolve().parent
    hashes={name:sha(root/name) for name in SOURCES}
    epoch=read(root/"START.json")["common_epoch"]
    a, ca, ga, controls=receipt_summary(local,hashes,epoch)
    b, cb, gb, others=receipt_summary(macbook,hashes,epoch)
    assert not controls and not others
    assert a["interval_half_open"][1] == b["interval_half_open"][0]
    fourway=intersect(intersect(ca,cb),intersect(ga,gb))
    pilot_a,pilot_b=[read(local/"pilot_{}.json".format(label)) for label in ("cpu","gpu")]
    assert all(pilot_a[key]==pilot_b[key] for key in KEYS)
    checksum=0
    witnesses=[]
    pilot_by_family=[Counter(), Counter()]
    start = pilot_a["first_id"]
    for i in range(start, start + pilot_a["count"]):
        row, cert = audit(i)
        if i in {start, start+2, start+8, start+10}:
            witnesses.append(cert)
        target = pilot_by_family[row[13]]
        target["cases"] += 1
        target["cases_with_actual_rank_loss"] += row[2] > row[3]
        target["actual_rank_loss_sum"] += row[2] - row[3]
        target["representative_collision_pairs"] += row[5]
        sig=mix(i)
        for v in row:
            sig=mix(sig^v)
        checksum^=sig
    assert f"{checksum:016x}"==pilot_a["positional_checksum64"]
    output.write_text(json.dumps(dict(common_epoch=epoch,local=a,macbook=b,
                                 independent_case_audits=pilot_a["count"], controls=witnesses,
                                 independent_pilot_diagnostics_by_family=dict(full=pilot_by_family[0], holed=pilot_by_family[1]),
                                 exact_expectations_all_linear_subspaces_and_cosets=[expected(start+8),expected(start+10)],
                                 four_way_active_seconds=sum(y-x for x,y in fourway),
                                 four_way_positive_intervals=len(fourway)),indent=2)+"\n")
    print(output)


if __name__=="__main__":
    main(*(Path(a) for a in sys.argv[1:]))
