"""Recompute the four-way overlap and in-window telemetry from compact input."""
import argparse
import bisect
import json
import statistics
from pathlib import Path

NS_PER_SECOND = 1_000_000_000
PERCENT_SCALE = 1_000_000


def union(intervals):
    result = []
    for start, end in sorted(intervals):
        assert start < end
        if result and start <= result[-1][1]:
            result[-1][1] = max(end, result[-1][1])
        else:
            result.append([start, end])
    return result


def intersection(left, right):
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        a,b = left[i]
        c,d = right[j]
        start,end = max(a,c), min(b,d)
        if start < end:
            result.append([start,end])
        if b < d:
            i += 1
        else:
            j += 1
    return result


def subset_samples(samples, ranges):
    starts = [pair[0] for pair in ranges]
    for sample in samples:
        k = bisect.bisect_right(starts,sample[0])-1
        if k >= 0 and sample[0] <= ranges[k][1]:
            yield sample


def stats(values, scale=1):
    values = list(values)
    assert values
    return {
        "samples": len(values), "minimum": min(values)/scale, "maximum": max(values)/scale,
        "mean": statistics.fmean(values)/scale, "median": statistics.median(values)/scale,
    }


def seconds(intervals):
    return sum(end-start for start,end in intervals)/NS_PER_SECOND


def verify(raw):
    assert raw["schema"] == "a28-r049-fourway-load-input-v1"
    math = raw["cpu_mathematical_run"]
    count_distribution = math["monochromatic_count_distribution"]
    assert sum(count_distribution.values()) == math["samples"]
    assert math["samples"] * math["affine_flats_per_sample"] == math["sample_flat_pairs_per_algorithm"]
    assert math["both_CPU_algorithms_agreed"] and math["first_three_incidence_moments_verified"]
    devices = raw["devices"]
    assert set(devices) == {"m6","m5-pro"}
    intervals = {}
    all_outer = []
    total_ANFs = 0
    total_coeffs = 0
    total_cpu_sampled_coeffs = 0
    identifier_ranges = []
    for name,d in devices.items():
        assert len(d["cpu_runs_ns_and_process_us"]) == 64
        assert d["logical_cpu_count"] == d["cpu_threads"]
        assert d["gpu_exit_code"] == d["cpu_exit_code"] == 0
        assert d["gpu_process_exited"] and d["cpu_processes_exited"]
        assert d["begin_candidate"] > 0 and d["candidate_count"] > 0
        assert d["candidate_count"]*512 == d["gpu_walsh_coefficients"]
        assert d["gpu_independently_sampled_coefficients"] % 512 == 0
        assert 0 < d["gpu_independently_sampled_coefficients"] <= d["gpu_walsh_coefficients"]
        cpu = union((start,end) for start,end,_ in d["cpu_runs_ns_and_process_us"])
        gpu = union(d["gpu_commands_union_ns"])
        assert gpu == d["gpu_commands_union_ns"]
        intervals[name] = (cpu,gpu)
        process = d["gpu_process_ns"]
        assert process[0] <= gpu[0][0] and gpu[-1][1] <= process[1]
        all_outer.extend(((cpu[0][0],cpu[-1][1]),process))
        identifier_ranges.append([d["begin_candidate"], d["begin_candidate"]+d["candidate_count"]])
        total_ANFs += d["candidate_count"]
        total_coeffs += d["gpu_walsh_coefficients"]
        total_cpu_sampled_coeffs += d["gpu_independently_sampled_coefficients"]
    assert max(x[0] for x in identifier_ranges) >= min(x[1] for x in identifier_ranges)
    combined_cpu = intersection(intervals["m6"][0],intervals["m5-pro"][0])
    combined_gpu = intersection(intervals["m6"][1],intervals["m5-pro"][1])
    four = intersection(combined_cpu,combined_gpu)
    outer = [[max(x[0] for x in all_outer),min(x[1] for x in all_outer)]]
    assert four and outer[0][0] < outer[0][1]
    assert all(outer[0][0] <= a <= b <= outer[0][1] for a,b in four)
    result = {}
    for name,d in devices.items():
        cpu_obs = list(subset_samples(d["cpu_sample_ns_active_pct_scaled_task_cores_scaled"],four))
        gpu_obs = list(subset_samples(d["gpu_sample_ns_whole_device_pct"],four))
        cpu_system = stats((obs[1] for obs in cpu_obs),PERCENT_SCALE)
        cpu_system["at_or_above_90_percent"] = sum(obs[1] >= 90*PERCENT_SCALE for obs in cpu_obs)
        task = stats((obs[2] for obs in cpu_obs),PERCENT_SCALE)
        gpu_system = stats(obs[1] for obs in gpu_obs)
        gpu_system["at_or_above_90_percent"] = sum(obs[1] >= 90 for obs in gpu_obs)
        result[name] = {
            "cpu_threads": d["cpu_threads"], "gpu_model": d["gpu_model"], "gpu_process_id": d["gpu_pid"],
            "own_cpu_core_equivalents_in_strict_intersection": task,
            "own_cpu_mean_fraction_of_logical_cores": task["mean"]/d["logical_cpu_count"],
            "whole_machine_cpu_samples_in_strict_intersection": cpu_system,
            "whole_device_gpu_samples_in_strict_intersection": gpu_system,
            "whole_sweep_cpu_interval_seconds": seconds(intervals[name][0]),
            "whole_sweep_own_metal_interval_seconds": seconds(intervals[name][1]),
        }
    return {
        "schema": "a28-r049-fourway-load-verified-v1",
        "nominal_outer_seconds": seconds(outer),
        "nominal_strict_active_seconds": seconds(four),
        "nominal_number_of_strict_active_intervals": len(four),
        "nominal_longest_strict_active_interval_seconds": max(b-a for a,b in four)/NS_PER_SECOND,
        "cpu_mathematical_samples": math["samples"],
        "cpu_sample_flat_pairs_per_algorithm": math["sample_flat_pairs_per_algorithm"],
        "cpu_samples_with_no_monochromatic_five_flat": count_distribution["0"],
        "total_gpu_distinct_ANFs": total_ANFs,
        "total_gpu_walsh_coefficients": total_coeffs,
        "total_exact_cpu_gpu_sampled_coefficients": total_cpu_sampled_coeffs,
        "devices": result,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=Path(__file__).with_name("input.json"))
    parser.add_argument("--output",type=Path)
    opts=parser.parse_args()
    computed = verify(json.loads(opts.input.read_bytes()))
    data = json.dumps(computed,indent=2,sort_keys=True)+"\n"
    if opts.output:
        opts.output.write_text(data)
    else:
        print(data,end="")


if __name__ == "__main__":
    main()
