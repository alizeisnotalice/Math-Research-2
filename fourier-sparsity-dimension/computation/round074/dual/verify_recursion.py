"""Independent exact Boolean recursion and certificate checking."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

SOURCES = ("recursion_cpu.c", "RecursionCommon.h", "RecursionKernel.metal", "RecursionMetal.swift", "run_host.py")
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
    assert diag[8]==sum(diag[23:39])==diag[16]+diag[17]==diag[9]+diag[10]
    assert diag[10]==diag[11]+diag[12]
    assert diag[15]==diag[13]-diag[14]
    csv_rows=list(csv.DictReader((path/"cpu_usage.csv").read_text().splitlines()))
    active=[r for r in csv_rows if int(r["active_workers"])==cpu["threads"]]
    machine_idle=[100-float(r["whole_machine_cpu_active_percent"]) for r in active
                  if float(r["whole_machine_cpu_active_percent"])>=0]
    cpu_actual=[float(r["a28_process_cpu_percent"])/100 for r in active]
    gpu_readings=[n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                  for n in sample["metrics"]["Device Utilization %"]]
    summary={key:cpu[key] for key in KEYS if key!="sums"}
    summary.update(
        one_bit_cases=2*cpu["count"]-diag[0], two_bit_cases=diag[0]-cpu["count"],
        total_parent_rank=diag[2], total_parent_points=diag[6], total_tail_points=diag[7],
        true_parent_bands=diag[8], repeated_bands=diag[9], exclusive_bands=diag[10],
        exclusive_inside_repeated_span=diag[11], exclusive_outside_repeated_span=diag[12],
        sum_child_tail_ranks=diag[13], total_true_tail_rank=diag[14], total_gluing_rank_defect=diag[15],
        uncertainty_equality_bands=diag[16], uncertainty_strict_bands=diag[17],
        cases_without_tail=diag[18], sum_of_missing_address_rank=diag[19],
        strict_quotient_gluing_cases=diag[22],
        band_histogram_by_active_children_and_parent_coefficients=[diag[j:j+4] for j in range(23,39,4)],
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


def walsh(values):
    if len(values)<=1:
        return values
    even=walsh(values[::2])
    odd=walsh(values[1::2])
    out=[]
    for a,b in zip(even,odd):
        out.extend([a+b,a-b])
    return out


def words(seed):
    state=seed or 1
    while True:
        state=(state^state<<13)&U32
        state^=state>>17
        state=(state^state<<5)&U32
        yield state


def parity(x):
    return x.bit_count()%2


def truth(case):
    k=1+(case>>2&1)
    n=5+(case>>3&1)
    family=case%4
    if case<6:
        k,n=(1 if case==0 else 2),5
    char_controls=((0,0),(4,4,1,5),(1,2,3,4),(0,1,2,4),(7,7,2,2),(1,1,1,2))
    descriptions=[]
    for b in range(1<<k):
        if case<6:
            shift=char_controls[case][b]
            phase=b%2 if case<4 else 0
            forms,table=[],0
            if case==0:
                forms,table=[1<<(j+2*b) for j in range(3)],1
            if case==3 and b==0:
                forms,table=[1,2],8
        else:
            rng=words(case^0xa2807401^(((b+1)*0x846ca68b)&U32))
            kind=next(rng)%4 if family==3 else family
            shift,phase=next(rng)%(1<<n),next(rng)%2
            target_rank=0 if kind==0 else 2 if kind==1 else 3+(next(rng)&1)
            table=0 if kind==0 else 8 if kind==1 else 1 if kind==2 else next(rng)%(1<<(1<<target_rank))
            forms=[]
            while len(forms)<target_rank:
                candidate=next(rng)%(1<<n)
                if rank(forms+[candidate])>len(forms):
                    forms.append(candidate)
        descriptions.append(dict(shift=shift,phase=phase,forms=forms,table=table))
    rows=[]
    for d in descriptions:
        phase,shift,table,forms=[d[key] for key in ("phase","shift","table","forms")]
        leaves=[]
        for x in range(1<<n):
            target=sum(parity(x&f)<<j for j,f in enumerate(forms))
            negative=(table>>target ^ phase ^ parity(shift&x))&1
            leaves.append(1-2*negative)
        rows.append(leaves)
    return k,n,family,rows,descriptions


def span(points):
    S={0}
    for x in points:
        S.update([v^x for v in S])
    return S


def audit_truth(k,n,family,truth_rows,certificate=False):
    B,Z=1<<k,1<<n
    child=[walsh(line) for line in truth_rows]
    parent=walsh([x for line in truth_rows for x in line])
    by_label={}
    for beta in range(Z):
        active=[b for b in range(B) if child[b][beta]]
        if active:
            by_label[beta]=active
    repeated=[beta for beta,t in by_label.items() if len(t)>1]
    core=span(repeated)
    rcore=rank(repeated)
    rdata=rank(by_label)
    assigned=[[] for _ in range(B)]
    for beta,sources in by_label.items():
        if beta not in core:
            assert len(sources)==1
            assigned[sources[0]].append(beta)
    per_rank=[rank(repeated+row)-rcore for row in assigned]
    rtail=rank(repeated+[beta for row in assigned for beta in row])-rcore
    assert rtail==rdata-rcore and rtail<=sum(per_rank)
    row=[0]*39
    row[0:2]=k,n
    row[5]=rcore
    row[8:10]=len(by_label),len(repeated)
    row[12]=sum(map(len,assigned))
    row[13:16]=sum(per_rank),rtail,sum(per_rank)-rtail
    row[21]=family
    row[22]=int(row[15]>0)
    row[18]=int(not row[12])
    supported,core_parent,glue_generators=[],[],[]
    audits=[]
    hash=0
    for beta in range(Z):
        active=by_label.get(beta,[])
        column=[parent[alpha*Z+beta] for alpha in range(B)]
        expected=[sum((-1 if parity(alpha&b) else 1)*child[b][beta] for b in active)
                  for alpha in range(B)]
        assert column==expected
        for x in column:
            hash=(hash*16777619&U32)^(x&U32)
        band=[alpha for alpha,c in enumerate(column) if c]
        size=len(band)
        if not active:
            assert not band
            continue
        a=len(active)
        assert a*size>=B
        row[16]+=a*size==B
        row[17]+=a*size>B
        row[23+4*(a-1)+size-1]+=1
        if a==1:
            assert size==B
            assert all(abs(val)==abs(child[active[0]][beta]) for val in column)
            row[10]+=1
            row[11]+=beta in core
        bits_here=[alpha*Z+beta for alpha in band]
        supported.extend(bits_here)
        if beta in core:
            core_parent.extend(bits_here)
        glue_generators.append(bits_here[0])
        glue_generators.extend(x^bits_here[0] for x in bits_here[1:])
        if certificate:
            audits.append(dict(beta=beta,active_children=active,unsigned_active_count=a,
                               child_coefficients=[x[beta] for x in child],full_parent_column=column,
                               parent_band_size=size,repeated_span=beta in core))
    d,dcore=rank(supported),rank(core_parent)
    assert d==rank(glue_generators)
    assert 0<=d-rdata<=k and 0<=dcore-rcore<=k
    assert not row[10] or d==rdata+k
    assert not row[11] or dcore==rcore+k
    assert not row[12] or d-dcore==rtail+k-(dcore-rcore)
    assert row[12] or d==dcore
    assert row[10]==row[11]+row[12]
    tail_points=len(supported)-len(core_parent)
    assert tail_points==B*row[12]
    row[2:5]=d,rdata,dcore
    row[6:8]=len(supported),tail_points
    row[19:21]=rdata+k-d,hash
    certificate_data=dict(k=k,leaf_dimensions=n,
                          truth_table_negative_masks=[sum((x<0)<<i for i,x in enumerate(rr)) for rr in truth_rows],
                          repeated_labels=repeated,span_repeated_labels=sorted(core),parent_rank=d,parent_core_rank=dcore,
                          child_tail_labels=assigned,individual_child_tail_quotient_ranks=per_rank,true_tail_quotient_rank=rtail,
                          signed_bands=audits,parent_tail_points=tail_points,diagnostics=row)
    return row,certificate_data if certificate else None


def audit_case(case,certificate=False):
    k,n,family,rows,desc=truth(case)
    row,cert=audit_truth(k,n,family,rows,certificate)
    if cert is not None:
        cert["id"]=case
        cert["leaf_descriptions"]=desc
        if k==2:
            top_specs=[]
            for first in range(2):
                order=[(a<<first)|(b<<(1-first)) for a in range(2) for b in range(2)]
                ordered_rows=[rows[v] for v in order]
                top,_=audit_truth(1,n+1,family,[ordered_rows[0]+ordered_rows[1],
                                               ordered_rows[2]+ordered_rows[3]])
                assert top[2]==row[2] and top[6]==row[6]
                children=[audit_truth(1,n,family,ordered_rows[j:j+2])[0] for j in (0,2)]
                top_specs.append(dict(first_address_bit=first,root_parent_rank=top[2],root_parent_size=top[6],
                                      intermediate_child_parent_ranks=[r[2] for r in children],
                                      intermediate_child_parent_sizes=[r[6] for r in children],
                                      root_tail_points=top[7],intermediate_tail_points=[r[7] for r in children]))
            cert["two_binary_recursion_orders"]=top_specs
    return row,cert


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
    total=[0]*39
    checksum=0
    examples=[]
    for case in range(p1["count"]):
        row,cert=audit_case(case,case<6)
        if cert:
            assert row==controls[str(case)]==p2["controls"][str(case)]
            examples.append(cert)
        sig=signature(case)
        for j,value in enumerate(row):
            total[j]+=value
            sig=signature(sig^value)
        checksum^=sig
    assert total==p1["sums"] and f"{checksum:016x}"==p1["positional_checksum64"]
    output.write_text(json.dumps(dict(common_epoch=epoch,local=a,macbook=b,
                      independent_pilot_cases=p1["count"],independent_pilot_bands=total[8],
                      explicit_controls=examples,four_way_overlap_seconds=sum(y-x for x,y in fourway),
                      four_way_positive_intervals=len(fourway)),indent=2)+"\n")
    print(output)


if __name__=="__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
