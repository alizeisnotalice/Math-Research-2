"""Independent Python generation, transforms, relations, and archive verification."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

SOURCES = ("clique_cpu.c", "AuditCommon.h", "ClickKernel.metal", "CliqueMetal.swift", "run_host.py")
KEYS = ("first_id", "count", "strict_graph", "strict_basis", "strict_activity",
        "sums", "positional_checksum64")
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


def relation_graph(labels, b):
    adj = [0] * (1 << b)
    for quad in combinations(labels, 4):
        a, y, z, t = quad
        if a ^ y ^ z ^ t == 0:
            clique = sum(1 << v for v in quad)
            for v in quad:
                adj[v] |= clique ^ (1 << v)
    return adj


def most_vertices(vertices, adj):
    """Bron-Kerbosch pivot search, with the clique itself as a witness."""
    best = []

    def descend(used, candidates, excluded):
        nonlocal best
        if len(used) + candidates.bit_count() <= len(best):
            return
        if not candidates:
            if not excluded and len(used) > len(best):
                best = used
            return
        pivot = max(bits(candidates | excluded), key=lambda v: (candidates & adj[v]).bit_count())
        for v in bits(candidates & ~adj[pivot]):
            descend(used + [v], candidates & adj[v], excluded & adj[v])
            candidates ^= 1 << v
            excluded |= 1 << v

    descend([], sum(1 << v for v in vertices), 0)
    return best


def tiny_clique_size(vertices, adj):
    if not vertices:
        return 0
    for size in range(len(vertices), 1, -1):
        for group in combinations(vertices, size):
            if all(adj[a] >> b & 1 for a, b in combinations(group, 2)):
                return size
    return 1


def compute(identifier, certificate=False, custom=None):
    a, b, family, tables = truth(identifier) if custom is None else custom
    B = 1 << b
    rows = [[1 - 2 * (table >> x & 1) for x in range(B)] for table in tables]
    supports = [tuple(x for x, v in enumerate(walsh(y)) if v) for y in rows]
    masks = [sum(1 << x for x in sup) for sup in supports]
    bands = [0] * B
    for t, v in enumerate(walsh([f for row in rows for f in row])):
        if v:
            bands[t % B] += 1
    labels = sorted(set().union(*supports))
    assert labels == [i for i, v in enumerate(bands) if v]
    q = sum(bands[1:])
    nonzero = [i for i in labels if i]
    r = rank(nonzero)
    actual = [0] * B
    for label_set in set(supports):
        local = relation_graph(label_set, b)
        for v, adj in enumerate(local):
            if v in label_set:
                assert set(bits(adj)) == set(label_set) - {v}
            actual[v] |= adj
    projected = relation_graph(sorted(set(labels) | {0}), b)
    assert all(a & ~p == 0 for a, p in zip(actual, projected))
    clique_a, clique_p = most_vertices(nonzero,actual), most_vertices(nonzero,projected)
    assert len(clique_a) <= len(clique_p)
    row = [0] * 37
    row[0], row[1], row[5], row[6] = r, q, len(clique_a), len(clique_p)
    row[7:10] = [r, r, r]
    row[10], row[12], row[13] = family, a, b
    packed = 0
    for mask, table in zip(masks, tables):
        packed = ((packed * 16777619) & U32) ^ (mask ^ table)
    for k in range(B):
        packed = ((packed * 16777619) & U32) ^ (bands[k] ^ (actual[k] << 8) ^ (projected[k] << 16))
    row[11] = packed
    row[14 + len(clique_p) - len(clique_a)] = 1
    best = {}
    best_activity = {}
    for I in combinations(nonzero, r):
        if rank(I) != r:
            continue
        pick = sum(1 << beta for beta in I)
        occ = max((pick & mask).bit_count() for mask in masks)
        wa, wp = tiny_clique_size(I,actual), tiny_clique_size(I,projected)
        budget = sum(bands[beta] for beta in I)
        assert r*r <= budget*occ <= q*occ and occ <= wa <= wp
        row[2] += 1
        row[3] = max(row[3], wa - occ)
        row[4] = max(row[4], wp - wa)
        row[7:10] = list(map(min,zip(row[7:10], [occ,wa,wp])))
        row[29 + wp - wa] += 1
        row[33 + wa - occ] += 1
        if certificate:
            key = (wp - wa, wa - occ, -occ)
            witness = dict(basis=list(I), occupancy=occ, actual=wa, projected=wp,
                           individual_basis_band_budget=budget)
            if not best or key > best["score"]:
                best = dict(witness, score=key)
            if not best_activity or wa-occ > best_activity["activity_gap"]:
                best_activity = dict(witness, activity_gap=wa-occ)
    assert row[2]
    if certificate:
        return row, dict(id=identifier, address_dimension=a, outer_dimension=b,
                         truth_table_masks=tables, slice_supports=supports, band_sizes=bands,
                         largest_actual_clique=clique_a, largest_projected_clique=clique_p,
                         designated_basis=best, largest_activity_gap_basis=best_activity, diagnostics=row)
    return row, None


def pair_plane_control():
    planes = list(combinations([1,2,4,8],2))
    tables = [sum(1 << x for x in range(16) if (u & x) and (v & x))
              for u,v in planes + planes[:2]]
    _, certificate=compute(-1,certificate=True,custom=(3,4,1,tables))
    basis=certificate["largest_activity_gap_basis"]
    assert basis["basis"]==[1,2,4,8] and basis["occupancy"]==2
    assert basis["actual"]==basis["projected"]==4
    return certificate


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
    assert sum(diag[14:29])==cpu["count"] and sum(diag[29:33])==sum(diag[33:37])==diag[2]
    assert sum(diag[15:29])==cpu["strict_graph"]
    csv_rows=list(csv.DictReader((path/"cpu_usage.csv").read_text().splitlines()))
    active=[r for r in csv_rows if int(r["active_workers"])==cpu["threads"]]
    machine_idle=[100-float(r["whole_machine_cpu_active_percent"]) for r in active
                  if float(r["whole_machine_cpu_active_percent"])>=0]
    cpu_actual=[float(r["a28_process_cpu_percent"])/100 for r in active]
    gpu_readings=[n for sample in receipt["GPU_samples"] if sample["gpu_running"]
                  for n in sample["metrics"]["Device Utilization %"]]
    summary={key:cpu[key] for key in KEYS if key!="sums"}
    summary.update(
        bases=diag[2], total_projected_rank=diag[0], total_q=diag[1],
        all_graph_clique_gap_histogram=diag[14:29],
        every_basis_clique_gap_histogram=diag[29:33],
        every_basis_actual_clique_minus_occupancy_histogram=diag[33:37],
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


def main(local, macbook, output):
    root=Path(__file__).resolve().parent
    inputs={key:sha(root/key) for key in SOURCES}
    epoch=read(root/"START.json")["common_epoch"]
    a,ca,ga,controls=host(local,inputs,epoch)
    b,cb,gb,others=host(macbook,inputs,epoch)
    assert a["interval_half_open"][0]==0 and a["interval_half_open"][1]==b["interval_half_open"][0]
    assert not others
    fourway=intersect(intersect(ca,cb),intersect(ga,gb))
    p1,p2=[read(local/f"pilot_{unit}.json") for unit in ("cpu","gpu")]
    assert all(p1[key]==p2[key] for key in KEYS)
    sum_row=[0]*37
    checksum,witnesses=0,[]
    for i in range(p1["count"]):
        row,cert=compute(i,certificate=(i<4))
        if cert:
            assert row==p2["controls"][str(i)]==controls[str(i)]
            witnesses.append(cert)
        sig=signature(i)
        for j,value in enumerate(row):
            sum_row[j]+=value
            sig=signature(sig^value)
        checksum^=sig
    assert sum_row==p1["sums"] and f"{checksum:016x}"==p1["positional_checksum64"]
    output.write_text(json.dumps(dict(common_epoch=epoch,local=a,macbook=b,
                                 independent_pilot_cases=p1["count"],independent_pilot_bases=sum_row[2],
                                 explicit_controls=witnesses, independent_pair_plane_control=pair_plane_control(),
                                 four_way_overlap_seconds=sum(y-x for x,y in fourway),
                                 four_way_positive_intervals=len(fourway)), indent=2)+"\n")
    print(output)


if __name__=="__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
