"""Verify the experiment and dependent sharp examples with exact Walsh arithmetic."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from functools import cache
from itertools import combinations
from pathlib import Path
from fractions import Fraction

SOURCES=("fiber_cpu.c","FiberCommon.h","FiberKernel.metal","FiberMetal.swift","run_host.py")
KEYS=("first_id","count","dependent_instances","sharp_instances","max_excess","max_fiber",
      "sums","positional_checksum64")
U32,U64=(1<<32)-1,(1<<64)-1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def sign(x):
    return 1-2*(x.bit_count()&1)


def bits(mask):
    while mask:
        place=mask&-mask
        yield place.bit_length()-1
        mask^=place


def rank(seq):
    pivots={}
    for v in seq:
        while v:
            leading=v.bit_length()
            if leading not in pivots:
                pivots[leading]=v
                break
            v^=pivots[leading]
    return len(pivots)


def wordstream(initial):
    state=initial or 1
    while True:
        state=(state^state<<13)&U32
        state^=state>>17
        state=(state^state<<5)&U32
        yield state


def model(case):
    v=2+(case>>2)%3
    y=1+(case>>4&1)
    if case==0:
        v,y=4,2
    rng=wordstream(case^0xa2806901)
    anchor=next(rng)%16
    directions=[]
    while len(directions)<4:
        d=next(rng)%16
        if rank(directions+[d])>len(directions):
            directions.append(d)
    return v,y,case%4,[anchor]+[anchor^d for d in directions]


def local(case,family,y,v):
    seed=case^0x691b156d^((y+1)*0x9e3779b9&U32)^((v+1)*0x846ca68b&U32)
    rng=wordstream(seed)
    omit=next(rng)%5
    coef=[0]*5
    if family==0:
        coef=[int(i!=omit) for i in range(5)]
    else:
        amount=2 if family==2 else 2*(next(rng)%3)
        for _ in range(amount):
            coef[next(rng)%5]+=1
    phases=next(rng)
    coef=[(-a if phases>>i&1 else a) for i,a in enumerate(coef)]
    return coef,bool(family>=2 and phases>>7&1),bool(phases>>8&1)


def transform(values):
    if len(values)<=1:
        return values
    even=transform(values[::2])
    odd=transform(values[1::2])
    out=[]
    for a,b in zip(even,odd):
        out.extend([a+b,a-b])
    return out


@cache
def private_truth(target,character,phase):
    if abs(target)==4:
        values=[target//4]*32
    elif abs(target)==2:
        values=[(target//2)*sign(int((z&3)==3)) for z in range(32)]
    else:
        values=[sign((z>>4) ^ (((z>>2)&1)*((z>>3)&1) if not character else 0) ^ phase)
                for z in range(32)]
    spec=transform(values)
    assert spec[0]%8==0 and spec[0]//8==target
    spectral=[i for i,v in enumerate(spec) if v]
    support=sum(1<<i for i in spectral)
    nonzero=[i for i in spectral if i]
    cost=Fraction(12*rank(nonzero)**2,len(nonzero)) if nonzero else Fraction(0)
    assert cost.denominator==1
    return spec[0]//8,support,int(cost)


@cache
def verify_bag(mask):
    labels=list(bits(mask|1))
    pairs=Counter(a^b for a,b in combinations(labels,2))
    for a,b,c,d in combinations(labels,4):
        if a^b^c^d==0:
            assert 2 in [pairs[a^b],pairs[a^c],pairs[a^d]]


def fingerprint(row,value):
    row[19]=(row[19]*16777619&U32)^(value&U32)


def audit(case, certificate=False):
    V,Y,family,labels=model(case)
    assert rank([a^labels[0] for a in labels[1:]])==4
    row=[0]*21
    row[1],row[2],row[20]=V,Y,family
    core=set()
    sharp_witness=None
    for y in range(1<<Y):
        g_all=[]
        support_y=0
        for v in range(1<<V):
            coef,char,phase=local(case,family,y,v)
            g,cost=[],[]
            for u in range(16):
                prescribed=sum(a*sign(k&u) for a,k in zip(coef,labels))
                actual,support,weight=private_truth(prescribed,char,phase)
                fingerprint(row,actual^support)
                support_y|=support
                g.append(actual)
                cost.append(weight)
            g_all.extend(g)
            projected={b for b,it in enumerate(transform(g)) if it}
            assert projected.issubset(labels)
            t=len(projected)
            assert t<=4
            upper=[{0:27,2:16,4:0}[abs(value)] for value in g]
            counts=Counter(abs(value) for value in g)
            assert all(a<=e for a,e in zip(cost,upper))
            e,a=192*t+sum(upper),192*t+sum(cost)
            assert a<=e<=1058
            if t==4:
                assert counts=={0:6,2:8,4:2} and e==1058
            row[0]+=1
            row[7]+=e==1058
            row[8]+=a==1058
            row[9+t]+=1
            row[14]=max(row[14],e)
            row[15]=max(row[15],a)
            row[16:19]=[row[16+i]+counts[2*i] for i in range(3)]
            if certificate and a==1058 and sharp_witness is None:
                sharp_witness=dict(y=y,v=v,g_times_four=g,active_projected_labels=sorted(projected),
                                   private_local_cost_times_twelve=cost,joint_cost=str(Fraction(a,192)))
        verify_bag(support_y)
        core.update(i for i,v in enumerate(transform(g_all)) if v)
    fibers=[[a>>4 for a in sorted(core) if a%16==b] for b in range(16)]
    for freqs in fibers:
        fingerprint(row,sum(1<<a for a in freqs))
    row[3],row[4]=len(core),rank([x^min(core) for x in core]) if core else 0
    row[5]=max(0,row[3]-1-row[4])
    row[6]=max(map(len,fibers))
    cert=dict(id=case,kernel_variables=V,address_variables=Y,projected_affine_labels=labels,
              common_support=sorted(core),fiber_sizes=list(map(len,fibers)),diagnostics=row,
              first_sharp_pair=sharp_witness) if certificate else None
    return row,cert


def large_control(vbits):
    labels=[1,2,4,8]
    g_all=[]
    minimum192=1058
    for v in range(1<<vbits):
        Q=sum((v>>i&1)*(v>>(i+1)&1) for i in range(0,vbits,2))
        coef=[sign((Q&1)^(v&i)) for i in range(4)]
        g,cost=[],[]
        bag=0
        for u in range(16):
            b=sum(c*sign(u&ell) for c,ell in zip(coef,labels))
            value,support,paid=private_truth(b,False,False)
            g.append(value)
            cost.append(paid)
            bag|=support
        verify_bag(bag)
        actual=[b for b,it in enumerate(transform(g)) if it]
        assert actual==labels
        assert Counter(map(abs,g))=={0:6,2:8,4:2}
        current=192*len(actual)+sum(cost)
        assert current==1058
        minimum192=min(minimum192,current)
        g_all.extend(g)
    core=[i for i,a in enumerate(transform(g_all)) if a]
    assert len(core)==4*(1<<vbits)
    assert sorted({x%16 for x in core})==labels
    r=rank([x^core[0] for x in core])
    assert r==vbits+3
    return dict(kernel_variables=vbits,projected_labels=labels,common_support_size=len(core),
                common_affine_rank=r,common_affine_dependence_excess=len(core)-1-r,
                tested_kernel_points=1<<vbits,constant_local_cost=str(Fraction(minimum192,192)))


def pair_hash(x):
    x=((x^(x>>30))*0xbf58476d1ce4e5b9)&U64
    x=((x^(x>>27))*0x94d049bb133111eb)&U64
    return x^(x>>31)


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
    assert diag[0]==sum(diag[9:14])
    assert sum(diag[16:19])==16*diag[0]
    assert diag[7]==diag[13] and diag[8]<=diag[7]
    csv_rows=list(csv.DictReader((path/"cpu_usage.csv").read_text().splitlines()))
    active=[r for r in csv_rows if int(r["active_workers"])==cpu["threads"]]
    machine_idle=[100-float(r["whole_machine_cpu_active_percent"]) for r in active
                  if float(r["whole_machine_cpu_active_percent"])>=0]
    cpu_actual=[float(r["a28_process_cpu_percent"])/100 for r in active]
    gpu_readings=[n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                  for n in sample["metrics"]["Device Utilization %"]]
    summary={key:cpu[key] for key in KEYS if key!="sums"}
    summary.update(
        address_kernel_pairs=diag[0], total_common_spectrum_size=diag[3],
        total_common_affine_rank=diag[4], case_coordinate_rank_excess=diag[5],
        local_shared_band_count_histogram=diag[9:14],
        exact_actual_sharp_pairs=diag[8], upper_cost_sharp_pairs=diag[7],
        zero_half_one_counts=diag[16:19],
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


def main(local_dir,macbook_dir,output):
    root=Path(__file__).resolve().parent
    hashes={name:sha(root/name) for name in SOURCES}
    epoch=read(root/"START.json")["common_epoch"]
    a,ca,ga,wa=host(local_dir,hashes,epoch)
    b,cb,gb,wb=host(macbook_dir,hashes,epoch)
    assert a["interval_half_open"][1]==b["interval_half_open"][0] and not wb
    pc,pg=[read(local_dir/f"pilot_{name}.json") for name in ("cpu","gpu")]
    assert all(pc[k]==pg[k] for k in KEYS)
    total=[0]*21
    signature=0
    examples=[]
    for i in range(pc["count"]):
        row,witness=audit(i,certificate=i<4)
        if witness:
            assert row==wa[str(i)]==pg["controls"][str(i)]
            examples.append(witness)
        code=pair_hash(i)
        for j,value in enumerate(row):
            total[j]+=value
            code=pair_hash(code^value)
        signature^=code
    assert total==pc["sums"] and f"{signature:016x}"==pc["positional_checksum64"]
    overlapping=intersect(intersect(ca,cb),intersect(ga,gb))
    data=dict(common_epoch=epoch,local=a,macbook=b,independent_pilot_cases=pc["count"],
              independent_pilot_pairs=total[0],explicit_controls=examples,
              independent_large_dependent_sharp_controls=[large_control(8),large_control(10)],
              four_way_overlap_seconds=sum(y-x for x,y in overlapping),
              four_way_positive_intervals=len(overlapping))
    output.write_text(json.dumps(data,indent=2)+"\n")
    print(output)


if __name__=="__main__":
    main(*(Path(x) for x in sys.argv[1:]))
