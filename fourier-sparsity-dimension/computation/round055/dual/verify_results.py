"""Reconcile CPU/Metal censuses, independently verify witnesses and summarize active loads."""
import argparse
import csv
import fractions
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def hash_file(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def walsh(table):
    return [sum((-1)**(((table>>x)&1)^((x&alpha).bit_count()%2)) for x in range(32))
            for alpha in range(32)]


def rank(vectors):
    pivot={}
    for a in vectors:
        while a:
            degree=a.bit_length()-1
            if degree in pivot:a^=pivot[degree]
            else:
                pivot[degree]=a
                break
    return len(pivot)


def independent_minimum(w,configs,target):
    bound=64-(64>>target)
    full=sum(abs(v) for v in w)
    support=[i for i,x in enumerate(w) if x]
    minimum=None
    witness=None
    for index,cfg in enumerate(configs):
        r=cfg["rank"]
        if r==target:
            for constant in range(1<<r):
                labels=[cfg["labels"][i]^constant for i in support]
                if rank(labels)!=target:continue
                gap=sum(abs(w[i]) for i,l in zip(support,labels) if l)-bound
                assert gap>=0
                if gap>0 and (minimum is None or gap<minimum):
                    witness={"projection_index":index,"projection_rank":r,"constant":constant,"fixed_phase_only_label":False}
                    minimum=gap
        if r==target-1:
            labels=[cfg["labels"][i]|(1<<r) for i in support]
            if rank(labels)!=target:continue
            gap=full-bound
            assert gap>=0
            if gap>0 and (minimum is None or gap<minimum):
                witness={"projection_index":index,"projection_rank":r,"fixed_phase_only_label":True}
                minimum=gap
    return minimum,witness


def union(intervals):
    result=[]
    for a,b in sorted(intervals):
        assert a<=b
        if result and a<=result[-1][1]:
            result[-1][1]=max(result[-1][1],b)
        else:result.append([a,b])
    return result


def intersection(a,b):
    i=j=0;result=[]
    while i<len(a) and j<len(b):
        l=max(a[i][0],b[j][0]);r=min(a[i][1],b[j][1])
        if l<r:result.append([l,r])
        if a[i][1]<=b[j][1]:i+=1
        else:j+=1
    return result


def seconds(intervals):
    return sum(b-a for a,b in intervals)


def validate_configs(path):
    cfg=json.loads(path.read_bytes())
    assert len(cfg)==373
    assert [sum(c["rank"]==r for c in cfg) for r in range(1,6)]==[31,155,155,31,1]
    kernels={}
    for c in cfg:
        labels=c["labels"]
        r=c["rank"]
        assert len(labels)==32 and set(labels)==set(range(1<<r))
        assert all(labels[x^y]==labels[x]^labels[y] for x in range(32) for y in range(32))
        kernel=tuple(i for i,z in enumerate(labels) if z==0)
        assert kernel not in kernels
        kernels[kernel]=r
    return cfg


def load(directory,cfg):
    receipt=json.loads((directory/"receipt.json").read_bytes())
    paths=[directory/"cpu.json",directory/"gpu.json"]
    c,g=map(lambda p:json.loads(p.read_bytes()),paths)
    assert receipt["remaining_task_process_pids"]==[] and receipt["cpu_gpu_exit_codes"]==[0,0]
    assert receipt["GPU_CPU_aggregate_exact_match"]
    assert receipt["cpu_gpu_outputs_sha256"]==[hash_file(p) for p in paths]
    assert receipt["sampled_cpu_telemetry_sha256"]==hash_file(directory/"cpu_usage.csv")
    for name,h in receipt["source_sha256"].items(): assert hash_file(ROOT/name)==h
    assert c["first_id"]==g["first_id"] and c["count"]==g["count"] and c["threads"]==receipt["logical_cores"]
    assert g["all_checks_passed"] and g["gpu_command_buffers"]
    assert len(c["worker_cpu_seconds"])==c["threads"] and min(c["worker_cpu_seconds"])>0
    ranks={}
    for r,row in c["by_rank"].items():
        raw=g["by_rank"][r]
        assert all(row[k]==v for k,v in raw.items())
        hist={int(k):v for k,v in row["gap_histogram_units32"].items()}
        assert sum(hist.values())+row["below_bound"]==row["full"]
        assert hist.get(0,0)==row["endpoint"]
        threshold=16 if int(r)>=4 else 8
        assert sum(v for k,v in hist.items() if 0<k<threshold)==row["positive_below_requested_gap"]
        gap=min(k for k in hist if k>0)
        assert gap==row["minimum_positive_units32"]
        independent,witness=independent_minimum(walsh(row["minimum_truth_mask"]),cfg,int(r))
        assert independent==gap
        ranks[r]={**raw,"minimum_positive_gap":str(fractions.Fraction(gap,32)),
                  "minimum_truth_mask":row["minimum_truth_mask"],"independent_minimum_witness":witness}
    samples=list(csv.DictReader((directory/"cpu_usage.csv").read_text().splitlines()))
    measured=[s for s in samples if int(s["active_workers"])==c["threads"]]
    assert measured
    occupied=union(g["gpu_command_intervals_epoch"])
    return {
        "device":g["device"],"cpu_pid":receipt["cpu_pid"],"metal_pid":receipt["metal_pid"],
        "start_epoch":receipt["same_requested_epoch"],"cpu_threads":c["threads"],
        "first_id":c["first_id"],"count":c["count"],"cpu_wall_seconds":c["wall"],
        "cpu_time_seconds":c["cpu_seconds"],"mean_a28_process_cpu_percent":100*c["cpu_seconds"]/c["wall"],
        "peak_a28_process_cpu_percent_in_all_active_thread_samples":receipt["peak_a28_cpu_process_percent"],
        "median_a28_process_cpu_percent_in_all_active_thread_samples":statistics.median(float(s["a28_process_cpu_percent"]) for s in measured),
        "whole_cpu_machine_peak_percent":receipt["peak_whole_machine_cpu_percent"],
        "a28_cpu_percent_samples_all_threads":len(measured),
        "sampled_whole_device_gpu_peak_percent":receipt["peak_whole_device_gpu_percent"],
        "sampled_device_gpu_while_metal_alive":sum(sum(len(v) for k,v in s["metrics"].items() if k=="Device Utilization %")
                                              for s in receipt["GPU_samples"] if s["gpu_running"]),
        "gpu_metal_command_seconds_sum":g["gpu_command_seconds_sum"],"gpu_metal_command_union_seconds":seconds(occupied),
        "gpu_integer_walsh_coefficients":g["count"]*32,
        "direct_CPU_Walsh_coefficients_crosschecked":g["direct_independent_cpu_coefficients"],
        "ranks":ranks,"cpu_gpu_exit_codes":[0,0],
    },[(c["start_epoch"],c["end_epoch"])],occupied


def strict_samples(directory, intervals):
    def inside(t): return any(lo <= t <= hi for lo,hi in intervals)
    cpu_rows = list(csv.DictReader((directory/"cpu_usage.csv").read_text().splitlines()))
    cpu = [x for x in cpu_rows if inside(float(x["unix_epoch"]))]
    related = [float(x["a28_process_cpu_percent"]) for x in cpu]
    machine = [float(x["whole_machine_cpu_active_percent"]) for x in cpu
               if float(x["whole_machine_cpu_active_percent"]) >= 0]
    receipt = json.loads((directory/"receipt.json").read_bytes())
    gpu = [n for x in receipt["GPU_samples"] if inside(x["epoch"])
           for n in x["metrics"]["Device Utilization %"]]
    def stats(values):
        return {"samples":len(values),"peak":max(values) if values else None,
                "mean":statistics.fmean(values) if values else None}
    return {"a28_process_cpu_percent":stats(related),
            "whole_machine_cpu_percent":stats(machine),"whole_machine_gpu_percent":stats(gpu)}


def main():
    parser=argparse.ArgumentParser()
    for name in ("local","macbook","out"):parser.add_argument(name,type=Path)
    a=parser.parse_args()
    cfg=validate_configs(ROOT/"configurations.json")
    m6,c0,g0=load(a.local,cfg)
    m5,c1,g1=load(a.macbook,cfg)
    assert m6["start_epoch"]==m5["start_epoch"]
    assert m6["first_id"]+m6["count"]<=m5["first_id"]
    simultaneous = intersection(intersection(c0,c1),intersection(g0,g1))
    assert simultaneous
    data={"schema":"a28-r055-independently-rechecked-four-way-census-v1",
          "strict_overlap_sampled_loads": {
              "local":strict_samples(a.local,simultaneous),
              "macbook":strict_samples(a.macbook,simultaneous),
          },
          "combined_total_full_pairs":sum(row["full"] for host in (m6,m5) for row in host["ranks"].values()),
          "combined_total_endpoints":sum(row["endpoint"] for host in (m6,m5) for row in host["ranks"].values()),
          "combined_gpu_integer_walsh_coefficients":m6["gpu_integer_walsh_coefficients"]+m5["gpu_integer_walsh_coefficients"],
          "combined_independent_walsh_coefficients":m6["direct_CPU_Walsh_coefficients_crosschecked"]+m5["direct_CPU_Walsh_coefficients_crosschecked"],
          "projections_by_rank":[31,155,155,31,1],"two_machine_truth_count":m6["count"]+m5["count"],
          "nominal_common_CPU_and_Metal_active_seconds":seconds(simultaneous),
          "nominal_two_CPU_overlap_seconds":seconds(intersection(c0,c1)),
          "nominal_two_Metal_overlap_seconds":seconds(intersection(g0,g1)),
          "CPU_and_GPU_percent_scopes_differ":True,
          "local":m6,"macbook":m5}
    a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    print(json.dumps(data,sort_keys=True))


if __name__=="__main__":main()
