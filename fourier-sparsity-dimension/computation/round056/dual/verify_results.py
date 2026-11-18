"""Audit the exact full five-variable scan and regenerate every structural certificate."""
import argparse
import collections
import csv
import hashlib
import json
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ANCHORS={0,1,2,4,8,16}
FREE=[x for x in range(32) if x not in ANCHORS]
AFFINE=[sum((((b&x).bit_count()^(b>>5))&1)<<x for x in range(32)) for b in range(64)]


def compact(mask):
    return sum(((mask>>x)&1)<<j for j,x in enumerate(FREE))


def expand(i):
    return sum(((i>>j)&1)<<x for j,x in enumerate(FREE))


def choose_char(mask):
    c=mask&1
    return (c<<5)|sum((((mask>>(1<<j))&1)^c)<<j for j in range(5))


def direct_walsh(mask):
    return [sum(1-2*(((mask>>x)&1)^((x&w).bit_count()&1)) for x in range(32)) for w in range(32)]


def affine_flat(m):
    if not m:return False
    p=[x for x in range(32) if (m>>x)&1]
    d={x^p[0] for x in p}
    return len(d).bit_count()==1 and all(x^y in d for x in d for y in d)


def flats_dictionary():
    spans={frozenset({0})}
    catalog={0:0}
    counts={}
    affine_counts={}
    for k in range(4):
        count=collections.Counter()
        codim=5-k
        for space in spans:
            cosets={frozenset(x^p for p in space) for x in range(32)}
            assert len(cosets)==1<<codim
            for coset in cosets:
                mask=sum(1<<x for x in coset)
                encoded=compact(mask^AFFINE[choose_char(mask)])
                if encoded in catalog: assert catalog[encoded]==codim
                catalog[encoded]=codim
                count[codim]+=1
        affine_counts[codim]=count[codim]
        counts[codim]=sum(value==codim for value in catalog.values())
        if k<3:
            spans={frozenset(space|{x^p for p in space}) for space in spans for x in range(32) if x not in space}
    assert affine_counts=={5:32,4:496,3:1240,2:620}
    assert counts=={5:32,4:496,3:1240,2:155} and len(catalog)==1924
    return catalog,counts,affine_counts


def union(intervals):
    out=[]
    for a,b in sorted(intervals):
        assert a<=b
        if out and a<=out[-1][1]:
            out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return out


def intersect(a,b):
    out=[]
    i=j=0
    while i<len(a) and j<len(b):
        x,y=a[i];u,v=b[j]
        if max(x,u)<min(y,v):out.append([max(x,u),min(y,v)])
        if y<=v:i+=1
        else:j+=1
    return out


def secs(intervals):
    return sum(y-x for x,y in intervals)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(folder):
    receipt=json.loads((folder/"receipt.json").read_bytes())
    source=[folder/"cpu.json",folder/"gpu.json"]
    cpu,gpu=map(lambda x:json.loads(x.read_bytes()),source)
    assert receipt["remaining_task_process_pids"]==[] and receipt["cpu_gpu_exit_codes"]==[0,0]
    assert receipt["GPU_CPU_aggregate_exact_match"]
    assert receipt["cpu_gpu_outputs_sha256"]==[sha(p) for p in source]
    assert receipt["sampled_cpu_telemetry_sha256"]==sha(folder/"cpu_usage.csv")
    for name,h in receipt["source_sha256"].items():assert sha(ROOT/name)==h
    for name in ("first_id","count","affine_variants_transformed_for_each_representative",
                 "distinct_raw_truth_tables_transformed","norm_histogram_units32","norm_position_checksum64"):
        assert cpu[name]==gpu[name]
    assert cpu["affine_variants_transformed_for_each_representative"]==64
    assert cpu["distinct_raw_truth_tables_transformed"]==64*cpu["count"]
    assert cpu["threads"]==receipt["logical_cores"] and min(cpu["worker_seconds"])>0
    assert sum(cpu["norm_histogram_units32"].values())==cpu["count"]
    assert gpu["walsh_coefficients_evaluated"]==32*cpu["distinct_raw_truth_tables_transformed"]
    assert gpu["all_checks_passed"] and [x[:2] for x in cpu["below_three_records"]]==gpu["below_three_id_and_norm"]
    assert sum(v for k,v in cpu["norm_histogram_units32"].items() if int(k)<96)==len(cpu["below_three_records"])
    sample=list(csv.DictReader((folder/"cpu_usage.csv").read_text().splitlines()))
    full=[float(v["a28_process_cpu_percent"]) for v in sample if int(v["active_workers"])==cpu["threads"]]
    assert full
    gpu_windows=union(gpu["gpu_command_intervals_epoch"])
    result={
        "device":gpu["device"],"cpus":cpu["threads"],"cpu_pid":receipt["cpu_pid"],"gpu_pid":receipt["metal_pid"],
        "input_half_open":[cpu["first_id"],cpu["first_id"]+cpu["count"]],"exit_codes":receipt["cpu_gpu_exit_codes"],
        "actual_raw_truth_tables_each_backend":cpu["distinct_raw_truth_tables_transformed"],
        "metal_walsh_coefficients":gpu["walsh_coefficients_evaluated"],
        "cpu_direct_crosschecked_coefficients":gpu["direct_CPU_walsh_coefficients_crosschecked"],
        "cpu_wall_seconds":cpu["wall"],"cpu_seconds":cpu["cpu_seconds"],
        "cpu_mean_process_percent":100*cpu["cpu_seconds"]/cpu["wall"],
        "cpu_sampled_peak_process_percent":max(full),
        "cpu_sampled_median_process_percent":statistics.median(full),
        "cpu_samples_all_workers":len(full),
        "whole_machine_cpu_peak_percent":receipt["peak_whole_machine_cpu_percent"],
        "whole_device_peak_gpu_percent":receipt["peak_whole_device_gpu_percent"],
        "metal_command_union_seconds":secs(gpu_windows),
        "gpu_wall_seconds":gpu["gpu_wall_seconds"],
        "start_epoch":receipt["same_requested_epoch"],
        "per_function_measured_norm_checksum":cpu["norm_position_checksum64"],
        "receipt_sha256":sha(folder/"receipt.json"),
    }
    return cpu,result,[(cpu["start_epoch"],cpu["end_epoch"])],gpu_windows


def certify(records,catalog):
    found={}
    witnesses={}
    for id,norm,codim,character,defect in records:
        assert id not in found
        assert id in catalog and catalog[id]==codim
        mask=expand(id)
        assert not any((mask>>x)&1 for x in ANCHORS)
        W=direct_walsh(mask)
        assert sum(abs(y) for y in W)==norm
        assert defect==mask^AFFINE[character]
        if codim==0: assert mask==defect==0 and norm==32
        else:
            assert 2<=codim<=5 and affine_flat(defect)
            assert defect.bit_count()==1<<(5-codim)
            assert norm==96-(128>>codim)
        found[id]=codim
        witnesses.setdefault(str(codim),
            {"id":id,"truth_mask_hex":f"{mask:08x}","norm_numerator_over_32":norm,
             "affine_character_index_low_bit_directions_high_bit_constant":character,
             "defect_mask_hex":f"{defect:08x}","defect_points":[i for i in range(32) if (defect>>i)&1],
             "true_walsh_integers":W,
             "distance_to_nearest_character":(32-max(abs(y) for y in W))//2})
    assert found==catalog
    return witnesses


def boundary_example(catalog):
    for id in range(1,8192):
        W=direct_walsh(expand(id))
        if sum(abs(w) for w in W)==96 and id not in catalog:
            mask=expand(id)
            assert all(not affine_flat(mask^c) for c in AFFINE)
            return {"id":id,"truth_mask_hex":f"{mask:08x}","norm_numerator_over_32":96,
                    "true_walsh_integers":W}
    raise ValueError("No boundary example found")


def main():
    p=argparse.ArgumentParser()
    for name in ("local","macbook","out"):p.add_argument(name,type=Path)
    args=p.parse_args()
    cat,codims,flats=flats_dictionary()
    a,x,ca,ga=read(args.local)
    b,y,cb,gb=read(args.macbook)
    assert x["start_epoch"]==y["start_epoch"]
    assert [x["input_half_open"],y["input_half_open"]]==[[0,25165824],[25165824,67108864]]
    counts=collections.Counter(a["norm_histogram_units32"])
    counts.update(b["norm_histogram_units32"])
    assert sum(counts.values())==1<<26
    witnesses=certify(a["below_three_records"]+b["below_three_records"],cat)
    assert len(a["below_three_records"])+len(b["below_three_records"])==1924
    spec={
        "schema":"a28-r056-independent-exhaustive-structure-summary-v1",
        "all_truth_tables_cpu":x["actual_raw_truth_tables_each_backend"]+y["actual_raw_truth_tables_each_backend"],
        "all_truth_tables_metal":x["actual_raw_truth_tables_each_backend"]+y["actual_raw_truth_tables_each_backend"],
        "metal_walsh_coefficients":x["metal_walsh_coefficients"]+y["metal_walsh_coefficients"],
        "independent_cpu_coefficient_checks":x["cpu_direct_crosschecked_coefficients"]+y["cpu_direct_crosschecked_coefficients"],
        "all_normalized_norm_histogram_units32":dict(sorted(counts.items(),key=lambda x:int(x[0]))),
        "below_three_patterns_normalized":{k:v for k,v in counts.items() if int(k)<96},
        "below_three_patterns_full_raw":{k:64*v for k,v in counts.items() if int(k)<96},
        "normalized_structural_counts_by_codimension":dict(sorted({0:1,**codims}.items())),
        "all_affine_flats_per_codimension_before_normalization":dict(sorted(flats.items())),
        "extremal_certificates":witnesses,"sharp_boundary_nonflat_certificate":boundary_example(cat),
        "nominal_two_cpu_overlap_seconds":secs(intersect(ca,cb)),
        "nominal_two_metal_overlap_seconds":secs(intersect(ga,gb)),
        "nominal_four_way_active_overlap_seconds":secs(intersect(intersect(ca,cb),intersect(ga,gb))),
        "local":x,"macbook":y,
    }
    args.out.write_text(json.dumps(spec,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:v for k,v in spec.items() if k not in ("extremal_certificates","all_normalized_norm_histogram_units32")},sort_keys=True))

if __name__=="__main__":main()
