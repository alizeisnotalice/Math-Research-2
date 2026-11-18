"""Fraction audit of weighted Walsh bands, and exact normal-allocation certificates."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from functools import cache
from fractions import Fraction
from itertools import combinations
from pathlib import Path

SOURCES = ("weight_cpu.c", "WeightCommon.h", "WeightKernel.metal", "WeightMetal.swift", "run_host.py")
KEYS = ("first_id", "count", "sums", "positional_checksum64")
U32, U64 = (1 << 32) - 1, (1 << 64) - 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def bits(mask):
    result = []
    while mask:
        low = mask & -mask
        result.append(low.bit_length() - 1)
        mask ^= low
    return result


def rank(seq):
    basis = {}
    for x in seq:
        while x:
            lead = x.bit_length()
            if lead not in basis:
                basis[lead] = x
                break
            x ^= basis[lead]
    return len(basis)


def signature(x):
    x = ((x ^ x >> 30) * 0xbf58476d1ce4e5b9) & U64
    x = ((x ^ x >> 27) * 0x94d049bb133111eb) & U64
    return x ^ x >> 31


def truth(identifier):
    a = 3 + (identifier >> 2 & 1)
    b = 3 + (identifier >> 3 & 1)
    family = identifier & 3
    if identifier == 0:
        a = b = 4
    elif identifier == 1:
        a = b = 3
    state = identifier ^ 0xa2806701 or 1

    def word():
        nonlocal state
        state = (state ^ state << 13) & U32
        state ^= state >> 17
        state = (state ^ state << 5) & U32
        return state

    B = 1 << b
    tables = []
    for y in range(1 << a):
        kind = word() % 3 if family == 3 else family
        shift, phase = word() % B, word() % 2
        u, v = 1, 2
        if kind == 1:
            u = word() % B
            while u == 0:
                u = word() % B
            v = word() % B
            while not v or v == u:
                v = word() % B
        if identifier == 0:
            kind, shift, phase = 0, y, 0
        if identifier == 1:
            kind, shift, phase = 1, 0, 0
            u, v = ((1, 2), (1, 4), (2, 4))[y % 3]
        if kind == 2:
            tables.append(word() % (1 << B))
            continue
        pattern = 0
        for x in range(B):
            value = ((shift & x).bit_count() ^ phase) & 1
            if kind == 1:
                value ^= ((u & x).bit_count() & 1) * ((v & x).bit_count() & 1)
            pattern |= value << x
        tables.append(pattern)
    return a, b, family, tables


def walsh(values):
    if len(values) == 1:
        return values
    evens, odds = walsh(values[::2]), walsh(values[1::2])
    result = []
    for a, b in zip(evens, odds):
        result.extend([a + b, a - b])
    return result


def concatenate(intervals):
    result = []
    for a,b in sorted(intervals):
        assert a <= b
        if result and a <= result[-1][1]:
            result[-1][1] = max(b,result[-1][1])
        else:
            result.append([a,b])
    return result


def intersect(a,b):
    i=j=0
    result=[]
    while i<len(a) and j<len(b):
        start,end=max(a[i][0],b[j][0]),min(a[i][1],b[j][1])
        if start<end:
            result.append([start,end])
        if a[i][1]<=b[j][1]:
            i+=1
        else:
            j+=1
    return result


def host(path, hashes, epoch):
    cpu, gpu, receipt = [read(path / x) for x in ("cpu.json", "gpu.json", "receipt.json")]
    assert all(cpu[key]==gpu[key] for key in KEYS)
    assert receipt["source_sha256"]==hashes and receipt["same_requested_epoch"]==epoch
    assert receipt["cpu_gpu_exit_codes"]==[0,0] and not receipt["remaining_task_process_pids"]
    assert receipt["cpu_gpu_outputs_sha256"]==[sha(path / x) for x in ("cpu.json", "gpu.json")]
    assert receipt["sampled_cpu_telemetry_sha256"]==sha(path / "cpu_usage.csv")
    start,end=receipt["input_interval_half_open"]
    assert cpu["first_id"]==start and cpu["count"]==end-start
    assert all(t>0 for t in cpu["worker_cpu_seconds"]) and cpu["threads"]==receipt["logical_cores"]
    assert gpu["all_checks_passed"] and gpu["gpu_command_buffers"]==len(gpu["gpu_command_intervals_epoch"])
    diag=cpu["sums"]
    for j in range(6,15,3):
        assert sum(diag[j:j+3])==cpu["count"]
    csv_rows=list(csv.DictReader((path/"cpu_usage.csv").read_text().splitlines()))
    active=[r for r in csv_rows if int(r["active_workers"])==cpu["threads"]]
    machine_idle=[100-float(r["whole_machine_cpu_active_percent"]) for r in active
                  if float(r["whole_machine_cpu_active_percent"])>=0]
    cpu_actual=[float(r["a28_process_cpu_percent"])/100 for r in active]
    gpu_readings=[n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                  for n in sample["metrics"]["Device Utilization %"]]
    summary={key:cpu[key] for key in KEYS if key!="sums"}
    summary.update(
        cases_with_holes=diag[5], sum_d=diag[0], sum_L=diag[1], sum_M=diag[2], sum_s=diag[3],
        products_vs_LM={label:diag[start:start+3] for label,start in
                        (("intrinsic",6),("varied_intrinsic",9),("inverse_activity",12))},
        strictly_improved_by_selected_cosets=diag[15:19],
        workers=cpu["threads"], cpu_wall_seconds=cpu["wall"], cpu_seconds=cpu["cpu_seconds"],
        project_cpu_average_equivalent_cores=cpu["cpu_seconds"]/cpu["wall"],
        project_cpu_median_equivalent_cores=statistics.median(cpu_actual) if cpu_actual else None,
        project_cpu_peak_equivalent_cores=max(cpu_actual, default=None),
        machine_cpu_idle_samples=len(machine_idle), machine_cpu_idle_zero_samples=machine_idle.count(0),
        machine_cpu_idle_median=statistics.median(machine_idle) if machine_idle else None,
        device=gpu["device"], gpu_buffers=gpu["gpu_command_buffers"],
        gpu_command_sum_seconds=gpu["gpu_command_seconds"], gpu_wall_seconds=gpu["metal_wall_seconds"],
        gpu_utilization_samples=len(gpu_readings),
        gpu_utilization_median=statistics.median(gpu_readings) if gpu_readings else None,
        gpu_utilization_peak=max(gpu_readings,default=None),
        gpu_baseline=receipt["GPU_before"]["metrics"]["Device Utilization %"],
        exit_codes=receipt["cpu_gpu_exit_codes"], remaining_pids=receipt["remaining_task_process_pids"],
        interval_half_open=[start,end])
    return summary, [[cpu["start_epoch"],cpu["end_epoch"]]], concatenate(gpu["gpu_command_intervals_epoch"]), gpu["controls"]


DEN=720720
RECIP=840
SCHEMES=("unweighted","intrinsic","varied_intrinsic","inverse_activity")


def spans(seq):
    result={0}
    for x in seq:
        result.update([s^x for s in result])
    return result


def seed_values(seed):
    current=seed or 1
    while True:
        current=(current^current<<13)&U32
        current^=current>>17
        current=(current^current<<5)&U32
        yield current


def signed_dot(a,b):
    return (a&b).bit_count()&1


def fold_hash(acc,value):
    return ((acc*16777619)&U32)^(value&U32)


def exact(identifier):
    a,b,family,tables=truth(identifier)
    A,B=1<<a,1<<b
    masks=[[1-2*(table>>z&1) for z in range(B)] for table in tables]
    local=[walsh(row) for row in masks]
    full=walsh([bit for row in masks for bit in row])
    d=rank([i for i,x in enumerate(full) if x])
    bands={}
    stream=seed_values(identifier^0xa2807201)
    for beta in range(B):
        live={t for t in range(A) if local[t][beta]}
        freq=[alpha for alpha in range(A) if full[alpha*B+beta]]
        if not freq:
            assert not live
            continue
        normal=spans([f^freq[0] for f in freq])
        origin=min(live)
        pack={t for t in range(A) if all(not signed_dot(t^origin,h) for h in normal)}
        assert pack.issubset(live) and len(pack)*len(normal)==A
        stretch=next(stream)%8+1
        intrinsic=Fraction(len(normal),len(freq))
        weights=[Fraction(1),intrinsic,intrinsic*stretch,Fraction(1,len(live))]
        bands[beta]=dict(actual=len(freq), filled=len(normal), freq=freq, direction=normal, live=live,
                         packing=pack, stretch=stretch, weights=weights)
    L=max(sum(t in band["live"] for band in bands.values()) for t in range(A))
    L_pack=max(sum(t in band["packing"] for band in bands.values()) for t in range(A))
    M=sum(band["filled"] for band in bands.values())
    S=sum(band["actual"] for band in bands.values())
    holes=any(band["actual"]<band["filled"] for band in bands.values())
    row=[d,L,M,S,L_pack,int(holes)]+[0]*14+[a,b,family]
    fingerprint=0
    for coef in full:
        fingerprint=fold_hash(fingerprint,coef)
    for band in bands.values():
        for item in (sum(1<<x for x in band["freq"]),sum(1<<x for x in band["packing"]),band["stretch"]):
            fingerprint=fold_hash(fingerprint,item)
    schemes=[]
    for i,key in enumerate(SCHEMES):
        overlap=max(sum(b["weights"][i] for b in bands.values() if t in b["live"]) for t in range(A))
        pack=max(sum(b["weights"][i] for b in bands.values() if t in b["packing"]) for t in range(A))
        reciprocal=sum(b["filled"]/b["weights"][i] for b in bands.values())
        assert d*d<=7*pack*reciprocal<=7*overlap*reciprocal
        quantity=overlap*reciprocal
        if i:
            compare=(quantity>L*M)-(quantity<L*M)+1
            row[6+3*(i-1)+compare]=1
        row[15+i]=int(pack<overlap)
        scales=[1,DEN,DEN,DEN]
        fact_scales=[1,1,RECIP,1]
        for item in (overlap*scales[i],pack*scales[i],reciprocal*fact_scales[i]):
            assert item.denominator==1
            fingerprint=fold_hash(fingerprint,int(item))
        schemes.append((overlap,pack,reciprocal,quantity))
    row[19]=fingerprint
    return row,schemes,bands,tables


@cache
def small_subspaces(a):
    initial=frozenset({0})
    known={initial}
    remaining=[initial]
    while remaining:
        old=remaining.pop()
        for x in range(1<<a):
            if x in old:
                continue
            extra=old|frozenset(u^x for u in old)
            if extra not in known:
                known.add(extra)
                remaining.append(extra)
    return sorted(known,key=lambda d:(len(d),sorted(d)))


@cache
def add_spaces(a,b):
    return frozenset(x^y for x in a for y in b)


def minimum_normal_cost(bands,dimension,scheme):
    initial=frozenset({0})
    target=frozenset(spans([x for band in bands.values() for x in band["direction"]]))
    states={initial:(Fraction(0),[])}
    for band in bands.values():
        allowable=[f for f in small_subspaces(dimension) if f.issubset(band["direction"])]
        cost_factor=band["weights"][scheme]/band["filled"]
        updated={}
        for used,(cost,choices) in states.items():
            for trial in allowable:
                if len(used&trial)!=1:
                    continue
                combined=add_spaces(used,trial)
                new_cost=cost+cost_factor*(len(trial).bit_length())**2
                if combined not in updated or new_cost<updated[combined][0]:
                    updated[combined]=(new_cost,choices+[sorted(trial)])
        states=updated
    return states[target]


def witness(identifier):
    row,weights,bands,tables=exact(identifier)
    prices={}
    for i,(overlap,pack,reciprocal,quantity) in enumerate(weights):
        optimum,assignment=minimum_normal_cost(bands,row[20],i)
        assert optimum<=7*pack
        prices[SCHEMES[i]]=dict(active_weighted_sup=str(overlap),coset_weighted_sup=str(pack),
                                reciprocal_sum=str(reciprocal),before_factor7=str(quantity),
                                minimum_normal_cost=str(optimum),allocated_normal_spaces=assignment)
    return dict(case=identifier,dimensions=row[20:22],d=row[0],L=row[1],M=row[2],s=row[3],
                true_boolean_row_masks=tables,LM=row[1]*row[2],packed_diagnostics=row,weights=prices,
                bands={beta:dict(n=data["actual"],m=data["filled"],active_addresses=sorted(data["live"]),
                                 selected_coset=sorted(data["packing"]),freq=data["freq"],stretch=data["stretch"])
                       for beta,data in bands.items()})


def main(local,macbook,output):
    root=Path(__file__).resolve().parent
    inputs={key:sha(root/key) for key in SOURCES}
    epoch=read(root/"START.json")["common_epoch"]
    a,ca,ga,controls=host(local,inputs,epoch)
    b,cb,gb,others=host(macbook,inputs,epoch)
    assert a["interval_half_open"][0]==0 and a["interval_half_open"][1]==b["interval_half_open"][0] and not others
    fourway=intersect(intersect(ca,cb),intersect(ga,gb))
    p1,p2=[read(local/f"pilot_{unit}.json") for unit in ("cpu","gpu")]
    assert all(p1[key]==p2[key] for key in KEYS)
    total=[0]*23
    checksum=0
    best,worst=(Fraction(1),0),(Fraction(1),0)
    for case in range(p1["count"]):
        row,schemes,_,_=exact(case)
        ratio=schemes[1][3]/schemes[0][3]
        if ratio<best[0]:
            best=(ratio,case)
        if ratio>worst[0]:
            worst=(ratio,case)
        sig=signature(case)
        for j,value in enumerate(row):
            total[j]+=value
            sig=signature(sig^value)
        checksum^=sig
        if str(case) in controls:
            assert row==controls[str(case)]==p2["controls"][str(case)]
    assert total==p1["sums"] and f"{checksum:016x}"==p1["positional_checksum64"]
    selected=[witness(x) for x in sorted({0,2,19,best[1],worst[1]})]
    output.write_text(json.dumps(dict(common_epoch=epoch,local=a,macbook=b,
                      independent_pilot_cases=p1["count"],
                      pilot_best_intrinsic_ratio=(str(best[0]),best[1]),pilot_worst_intrinsic_ratio=(str(worst[0]),worst[1]),
                      witnesses=selected,four_way_overlap_seconds=sum(y-x for x,y in fourway),
                      four_way_positive_intervals=len(fourway)),indent=2)+"\n")
    print(output)


if __name__=="__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
