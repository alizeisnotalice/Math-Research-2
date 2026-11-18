"""Exact affine-flat geometry and signed-cut audit for the R32/R33 example."""
import collections
import ctypes
import datetime
import hashlib
import itertools
import json
import multiprocessing
import os
import pathlib
import platform
import resource
import subprocess
import sys
import threading
import time

SIZE = 128
BITCOUNT = [bin(i).count("1") for i in range(SIZE)]
POINT_BITS = [1 << i for i in range(SIZE)]


def span(basis):
    points = [0]
    for row in basis:
        points.extend([point ^ row for point in points])
    return points


def rref_chunks():
    """Each direction has exactly one RREF basis, with low bit in column zero."""
    tasks = []
    for dimension in range(8):
        for pivots in itertools.combinations(range(7), dimension):
            indices = [i for i in range(7) if i not in pivots]
            free = [(row, column) for row, pivot in enumerate(pivots)
                    for column in indices if column > pivot]
            for first in range(0, 1 << len(free), 32):
                tasks.append((dimension, pivots, indices, free, first,
                              min(first + 32, 1 << len(free))))
    return tasks


def local_walsh(values):
    result = values[:]
    stride = 1
    while stride < len(result):
        for start in range(0, len(result), 2 * stride):
            for index in range(start, start + stride):
                left, right = result[index], result[index + stride]
                result[index], result[index + stride] = left + right, left - right
        stride *= 2
    return result


def initialize(points, negative_mask, dual_walsh):
    global TABLES
    physical = [0] * SIZE
    for index, point in enumerate(points):
        physical[point] = -1 if negative_mask >> index & 1 else 1
    # On the ambient eight-dimensional dual space, this sign extends
    # identically across the free high bit. The Walsh cutoff on F2^7 is 8.
    dual = [value // 16 for value in dual_walsh]
    TABLES = ((physical, 16), (dual, 8))


def count_chunk(task):
    dimension, pivots, indices, free, first, stop = task
    started = time.process_time()
    histograms = [collections.Counter(), collections.Counter()]
    cut_histograms = [collections.Counter(), collections.Counter()]
    contained = [[], []]
    legal_proper_cuts = [set(), set()]
    maximum = [0, 0]
    maximum_multiplicity = [0, 0]
    main_sizes = [sum(value != 0 for value in table) for table, _ in TABLES]
    translates = span([1 << i for i in indices])
    ground_basis = [1 << pivot for pivot in pivots]
    for word in range(first, stop):
        basis = ground_basis[:]
        for bit, (row, column) in enumerate(free):
            if word >> bit & 1:
                basis[row] ^= 1 << column
        direction = span(basis)
        for translate in translates:
            points = [translate ^ point for point in direction]
            for index, (signs, cutoff) in enumerate(TABLES):
                values = [signs[point] for point in points]
                weight = sum(value != 0 for value in values)
                histograms[index][weight] += 1
                if weight > maximum[index]:
                    maximum[index], maximum_multiplicity[index] = weight, 0
                if weight == maximum[index]:
                    maximum_multiplicity[index] += 1
                if weight == len(points) and dimension >= 2:
                    contained[index].append((tuple(basis), translate))
                if not weight:
                    continue
                spectrum = local_walsh(values)
                spectral_magnitude_histogram = tuple(sorted(collections.Counter(map(abs, spectrum)).items()))
                cut_histograms[index][(weight, spectral_magnitude_histogram)] += 1
                if weight < main_sizes[index] and all(abs(value) in (0, cutoff) for value in spectrum):
                    legal_proper_cuts[index].add(sum(POINT_BITS[p] for p, value in zip(points, values) if value))
    return (dimension, stop - first, histograms, cut_histograms, contained,
            legal_proper_cuts, maximum, maximum_multiplicity, os.getpid(), time.process_time() - started)


def cpu_snapshot():
    """Host-wide macOS CPU ticks measured via HOST_CPU_LOAD_INFO."""
    libc = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
    libc.mach_host_self.restype = ctypes.c_uint
    values = (ctypes.c_uint * 4)()
    count = ctypes.c_uint(4)
    rc = libc.host_statistics(libc.mach_host_self(), 3, ctypes.byref(values), ctypes.byref(count))
    if rc != 0 or count.value != 4:
        raise RuntimeError("host_statistics failed: " + str(rc))
    return list(values)


def monitor_cpu(stop, samples):
    previous = cpu_snapshot()
    previous_time = time.monotonic()
    while not stop.wait(0.1):
        current = cpu_snapshot()
        stamp = time.monotonic()
        delta = [new - old for new, old in zip(current, previous)]
        if sum(delta) > 0:
            samples.append({
                "time_monotonic": stamp, "interval_seconds": stamp - previous_time,
                "host_idle_percent": round(100 * delta[2] / sum(delta), 2),
                "ticks": delta,
            })
        previous, previous_time = current, stamp


def geometry_record(chunks, index, dimension):
    intersections, signed_spectra = collections.Counter(), collections.Counter()
    max_weight, how_many = 0, 0
    full_flats, legal_cuts = [], set()
    for chunk in chunks:
        _, _, hist, cut_hist, contained, cuts, peak, multiplicity, _, _ = chunk
        intersections.update(hist[index])
        signed_spectra.update(cut_hist[index])
        full_flats.extend(contained[index])
        legal_cuts.update(cuts[index])
        if peak[index] > max_weight:
            max_weight, how_many = peak[index], 0
        if peak[index] == max_weight:
            how_many += multiplicity[index]
    result = {
        "intersection_histogram": dict(sorted(intersections.items())),
        "max_intersection": max_weight,
        "number_attaining_max_intersection": how_many,
        "contained_flats_count": intersections[1 << dimension],
        "contained_flats": [{"basis": basis, "translate": translate}
                            for basis, translate in sorted(full_flats)],
        "distinct_walsh_compatible_proper_cuts": [hex(cut) for cut in sorted(legal_cuts)],
        "nonempty_cut_local_walsh_absolute_histogram": [
            {"cut_weight": weight, "counts_by_abs_value": dict(walsh_histogram), "multiplicity": count}
            for (weight, walsh_histogram), count in sorted(signed_spectra.items())
        ],
    }
    return result


def gaussian_binomial(n, k):
    numerator, denominator = 1, 1
    for i in range(k):
        numerator *= (1 << n) - (1 << i)
        denominator *= (1 << k) - (1 << i)
    assert numerator % denominator == 0
    return numerator // denominator


def run(input_file, metal_file, output_file, workers):
    input_path, metal_path, output_path = map(pathlib.Path, (input_file, metal_file, output_file))
    signed_input = json.loads(input_path.read_text())
    points = signed_input["points"]
    negative_mask = signed_input["solutions"][0]
    dual_walsh = [
        sum((1 - 2 * ((negative_mask >> index) & 1)) * (1 - 2 * (BITCOUNT[frequency & point] & 1))
            for index, point in enumerate(points))
        for frequency in range(SIZE)
    ]
    if dual_walsh != json.loads(metal_path.read_text())["reference_walsh"]:
        raise ValueError("Independent core Walsh spectrum disagrees with Metal result")
    tasks = rref_chunks()
    # Large-span cases first, so tail effects are limited.
    tasks.sort(key=lambda task: -(task[5] - task[4]))
    ctx = multiprocessing.get_context("fork")
    samples, stop = [], threading.Event()
    with ctx.Pool(workers, initializer=initialize, initargs=(points, negative_mask, dual_walsh)) as pool:
        pids = [process.pid for process in pool._pool]
        monitor = threading.Thread(target=monitor_cpu, args=(stop, samples), daemon=True)
        monitor.start()
        start_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()
        start = time.monotonic()
        try:
            audit = list(pool.imap_unordered(count_chunk, tasks, chunksize=1))
        finally:
            end = time.monotonic()
            stop.set()
            monitor.join()
        pool.close()
        pool.join()

    dimensions = {}
    children_pid_lines = subprocess.check_output(["/bin/ps", "-Ao", "pid="], text=True)
    live = set(map(int, children_pid_lines.split()))
    observed = [sample for sample in samples if sample["time_monotonic"] <= end]
    for dimension in range(8):
        these_chunks = [chunk for chunk in audit if chunk[0] == dimension]
        directions = sum(chunk[1] for chunk in these_chunks)
        assert directions == gaussian_binomial(7, dimension)
        cosets = directions * (SIZE >> dimension)
        physical, dual = (geometry_record(these_chunks, index, dimension) for index in range(2))
        for record, cardinality in ((physical, 32), (dual, 16)):
            histogram = record["intersection_histogram"]
            assert sum(histogram.values()) == cosets
            assert sum(i * value for i, value in histogram.items()) == directions * cardinality
            assert sum(i * (i - 1) * value for i, value in histogram.items()) == (
                cardinality * (cardinality - 1) * (gaussian_binomial(6, dimension - 1) if dimension else 0)
            )
        dimensions[dimension] = {
            "directions": directions, "affine_flats": cosets, "S32": physical, "D16": dual,
        }
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    worker_seconds = collections.Counter()
    for chunk in audit:
        worker_seconds[chunk[-2]] += chunk[-1]
    result = {
        "schema": "a28-r034-affine-flat-audit-v1",
        "signed_input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "metal_input_sha256": hashlib.sha256(metal_path.read_bytes()).hexdigest(),
        "program_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "start_utc": start_utc,
        "platform": platform.platform(),
        "python_version": sys.version,
        "workers": workers,
        "logical_cores": os.cpu_count(),
        "rref_chunks": len(tasks),
        "worker_pids": pids,
        "worker_exit_verified": not (set(pids) & live),
        "active_workers": len(worker_seconds),
        "worker_task_cpu_seconds_by_pid": dict(worker_seconds),
        "worker_task_cpu_seconds_total": sum(worker_seconds.values()),
        "worker_cpu_seconds_including_setup_and_teardown": usage.ru_utime + usage.ru_stime,
        "parallel_wall_seconds": end - start,
        "parallel_average_core_equivalents": sum(worker_seconds.values()) / (end - start),
        "host_cpu_idle_samples": observed,
        "min_host_idle_pct": min((sample["host_idle_percent"] for sample in observed), default=None),
        "dimensions": dimensions,
    }
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    summary = {k: v for k, v in result.items() if k not in ("dimensions", "host_cpu_idle_samples")}
    for dimension in (2, 3, 4):
        summary[dimension] = {
            name: dimensions[dimension][name]["intersection_histogram"] for name in ("S32", "D16")
        }
    print(json.dumps(summary, sort_keys=False))


if __name__ == "__main__":
    run(*sys.argv[1:4], workers=int(sys.argv[4]))
