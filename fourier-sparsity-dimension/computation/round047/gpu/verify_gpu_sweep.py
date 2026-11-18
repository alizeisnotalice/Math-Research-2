"""Audit sweep histograms, candidate ranges, timing, and a few CPU sample ANFs."""
import argparse
import collections
import hashlib
import json
import pathlib
import struct
import subprocess

COEFFS = 512
WORD_MASK = (1 << 32) - 1


def mix(x):
    x ^= x >> 16
    x = x * 0x7FEB352D & WORD_MASK
    x ^= x >> 15
    x = x * 0x846CA68B & WORD_MASK
    return x ^ (x >> 16)


def cpu_spectrum(sample, seed):
    monomials = [x for x in range(COEFFS) if x.bit_count() == 3]
    monomials += [x for x in range(COEFFS) if x.bit_count() < 3]
    words = [sample] + [mix(sample ^ seed ^ salt) for salt in
                        (0x9E3779B9, 0x3C6EF372, 0xDAA66D2B, 0x78DDE6E4)]
    # Evaluate only the ANF coefficients for monomials present at the given x.
    terms = [term for i, term in enumerate(monomials) if (words[i//32] >> (i % 32)) & 1]
    values = [(1 - 2 * (sum(x & term == term for term in terms) % 2)) for x in range(COEFFS)]
    for step in (1,2,4,8,16,32,64,128,256):
        for start in range(0,COEFFS,2*step):
            for x in range(start,start+step):
                a,b = values[x],values[x+step]
                values[x],values[x+step] = a+b,a-b
    return values


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_union(intervals):
    result = []
    for a,b in sorted(intervals):
        assert a <= b
        if result and a <= result[-1][1]:
            result[-1] = (result[-1][0], max(b,result[-1][1]))
        else:
            result.append((a,b))
    return result


def duration(intervals):
    return sum(b-a for a,b in intervals)


def intersections(left, right):
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        a,b = left[i]
        c,d = right[j]
        if min(b,d) > max(a,c):
            result.append((max(a,c),min(b,d)))
        if b < d:
            i += 1
        else:
            j += 1
    return result


def check_device(directory, sample_gpu_binary):
    result_path = directory / "cubic_sweep_result.json"
    receipt_path = directory / "cubic_sweep_receipt.json"
    result = json.loads(result_path.read_bytes())
    receipt = json.loads(receipt_path.read_bytes())
    assert receipt["exit_code"] == 0 and receipt["process_exited"]
    assert receipt["result_sha256"] == digest(result_path)
    assert receipt["source"]["source_sha256"] == digest(pathlib.Path(__file__).parent / "R47CubicWalshLoad.swift")
    assert receipt["runner_sha256"] == digest(pathlib.Path(__file__).parent / "monitor_cubic.py")
    assert result["pid"] == receipt["process_id"]
    assert result["parseval_and_sample_checks_passed"]
    n = result["candidate_count"]
    first = result["first_candidate_id"]
    assert result["walsh_coefficients_evaluated"] == n*COEFFS
    absolute = result["absolute_walsh_histogram_even_bins_0_to_512"]
    sparse = result["sparsity_histogram"]
    assert len(absolute) == 257 and len(sparse) == 513
    assert sum(absolute) == n*COEFFS and sum(sparse) == n
    assert sum(count * (2*i)**2 for i,count in enumerate(absolute)) == n*COEFFS**2
    assert sum(s*count for s,count in enumerate(sparse)) == n*COEFFS - absolute[0]
    sampled = result["sampled_ids"]
    assert len(sampled) == len(set(sampled)) and all(first <= i < first+n for i in sampled)
    assert result["sampled_exact_cpu_gpu_coefficients"] == len(sampled)*COEFFS
    seed = int(result["seed_hex"],16)
    rows = []
    # Record reproducible Python reference samples separate from the on-device Swift direct sums.
    for sample in sorted(sampled)[::max(1,len(sampled)//8)][:8]:
        spectrum = cpu_spectrum(sample,seed)
        abs_counts = collections.Counter(abs(c) for c in spectrum)
        assert sum(c*c for c in spectrum) == COEFFS**2
        assert sum(c != 0 for c in spectrum) > 0
        entry = {
            "id": sample, "sparsity": sum(c != 0 for c in spectrum),
            "coefficient_sha256_le_int16": hashlib.sha256(struct.pack("<512h",*spectrum)).hexdigest(),
            "absolute_histogram": dict(sorted(abs_counts.items())),
        }
        if sample_gpu_binary is not None and len(rows) < 2:
            sample_file = directory / f"single_gpu_witness_{sample}.json"
            subprocess.run([str(sample_gpu_binary.resolve()),str(sample_file.resolve()),str(sample),
                            "1",result["seed_hex"],"256","1"], check=True)
            single = json.loads(sample_file.read_bytes())
            assert single["sampled_exact_cpu_gpu_coefficients"] == COEFFS
            assert single["absolute_walsh_histogram_even_bins_0_to_512"] == [
                abs_counts[2*i] for i in range(257)
            ]
            entry["single_gpu_histogram_crosscheck_sha256"] = digest(sample_file)
        rows.append(entry)
    intervals = interval_union(result["gpu_command_intervals_utc_unix"])
    data = [item["measurements"]["Device Utilization %"] for item in receipt["sampled_while_task_alive"]]
    recent = [x for row in data for x in row]
    assert sum(n>=90 for n in recent) == receipt["system_samples_at_least_90_percent"]
    return {
        "device": result["device"], "pid": result["pid"],
        "candidate_interval": [first,first+n], "number_of_anfs": n,
        "walsh_coefficients": n*COEFFS,
        "sparsity_min": min(i for i,v in enumerate(sparse) if v),
        "sparsity_max": max(i for i,v in enumerate(sparse) if v),
        "exact_sampled_coefficients_in_swift": result["sampled_exact_cpu_gpu_coefficients"],
        "time_interval": [receipt["observed_start_epoch"],receipt["observed_end_epoch"]],
        "own_metal_command_union_seconds": duration(intervals),
        "sum_individual_command_seconds": result["gpu_command_buffer_seconds_total"],
        "number_of_buffers": result["gpu_command_buffers"],
        "system_utilization_samples": len(recent),
        "system_samples_at_or_above_90pct": sum(x>=90 for x in recent),
        "system_samples_at_100pct": sum(x==100 for x in recent),
        "system_device_range_pct": [min(recent),max(recent)] if recent else None,
        "baseline_system_device_pct": receipt["baseline_ioaccelerator"]["measurements"]["Device Utilization %"],
        "exit_code": receipt["exit_code"], "process_exited": receipt["process_exited"],
        "independent_python_cpu_references": rows,
        "result_sha256": digest(result_path), "receipt_sha256": digest(receipt_path),
    }, intervals


def main():
    p = argparse.ArgumentParser()
    p.add_argument("local",type=pathlib.Path)
    p.add_argument("macbook",type=pathlib.Path)
    p.add_argument("output",type=pathlib.Path)
    p.add_argument("--sample-gpu-binary",type=pathlib.Path)
    args = p.parse_args()
    a,a_gpu = check_device(args.local,args.sample_gpu_binary)
    b,b_gpu = check_device(args.macbook,args.sample_gpu_binary)
    assert a["candidate_interval"][1] <= b["candidate_interval"][0] or b["candidate_interval"][1] <= a["candidate_interval"][0]
    overlap = intersections([tuple(a["time_interval"])],[tuple(b["time_interval"])])
    raw = {
        "schema": "a28-r047-sustained-dual-gpu-verified-summary-v1",
        "local": a, "macbook": b,
        "total_anfs": a["number_of_anfs"]+b["number_of_anfs"],
        "total_walsh_coefficients": a["walsh_coefficients"]+b["walsh_coefficients"],
        "nominal_overlapping_process_seconds": duration(overlap),
        "nominal_overlapping_metal_command_seconds": duration(intersections(a_gpu,b_gpu)),
        "system_telemetry_includes_other_processes": True,
    }
    args.output.write_text(json.dumps(raw,indent=2,sort_keys=True)+"\n")
    print({k:v for k,v in raw.items() if k not in ("local","macbook")})


if __name__ == "__main__":
    main()
