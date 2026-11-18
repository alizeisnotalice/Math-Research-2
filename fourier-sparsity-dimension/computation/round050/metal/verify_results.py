"""Check two receipts, finite corpus coverage, and independent direct Walsh witnesses."""
import argparse
import hashlib
import json
import struct
from pathlib import Path

LINES=((0,1),(2,4),(3,7),(5,6))


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def endpoint_specs(line,choice):
    label=line ^ ((choice%4)<<2)
    a=(1,1,4,2)[line]
    b=(2,4,8,8)[line]
    o=(4,2,1,1)[line]
    c=lambda l:(l,0,0,0,0)
    q=lambda l,r,s,sign=0,shift=0:(l,r,s,sign,shift)
    modes=(
        (c(label),c(label)),
        (c(label),c(label^o)),
        (q(label,a,b),q(label^a,a,b,1)),
        (q(label,a,b),q(label^o,a,b)),
        (q(label,a,b),q(label,a^o,b,0,1)),
        (c(label),q(label,a,b)),
        (q(label,a,b),c(label^o)),
        (q(label^b,a^o,b),q(label^a,a,b,0,1)),
    )
    return modes[choice//4]


def value(y,spec):
    label,a,b,constant,shift=spec
    dot=lambda x,z:(x&z).bit_count()%2
    return (-1)**(dot(y^shift,a)*dot(y^shift,b)^dot(y,label)^constant)


def table(case_id):
    result=[0]*128
    for i,address in enumerate(LINES):
        inputs=endpoint_specs(i,(case_id>>(5*i))&31)
        for t in range(2):
            for y in range(16):
                result[address[t]+8*y]=value(y,inputs[t])
    return result


def walsh(truth):
    return [sum((1-2*((x&alpha).bit_count()%2))*value for x,value in enumerate(truth))
            for alpha in range(len(truth))]


def verify_witness(witness):
    f=table(witness["case_id"])
    beta=witness["data_frequency"]
    D=[[sum((1-2*((y&beta).bit_count()%2))*f[a+8*y] for y in range(16)) for a in line]
       for line in LINES]
    local=[(a+b!=0)+(a-b!=0) for a,b in D]
    got=walsh(f)[8*beta:8*beta+8]
    assert got == witness["address_walsh"]
    assert local == witness["local_preimages"]
    assert sum(v!=0 for v in local)==witness["active_lines"]
    assert sum(int(v!=0) for pair in D for v in pair)==witness["supported_addresses"]
    assert 4*sum(local)==witness["precollision_occurrences"]
    assert sum(v!=0 for v in got)==witness["true_band_sparsity"]
    union=0
    for alpha in range(8):
        if any(sum((1-2*((alpha&a).bit_count()%2))*v for a,v in zip(line,d)) != 0
               for line,d in zip(LINES,D)):
            union+=1
    assert union==witness["precollision_union"]
    return witness["case_id"]


def merge(intervals):
    result=[]
    for a,b in sorted(intervals):
        assert a<=b
        if result and a<=result[-1][1]:
            result[-1][1]=max(result[-1][1],b)
        else:
            result.append([a,b])
    return result


def intersection_time(a,b):
    i=j=0
    total=0.0
    while i<len(a) and j<len(b):
        x,y=a[i]
        u,v=b[j]
        total+=max(0,min(y,v)-max(x,u))
        if y<=v: i+=1
        else: j+=1
    return total


def read(directory):
    result_file=directory/"result.json"
    data=result_file.read_bytes()
    result=json.loads(data)
    receipt=json.loads((directory/"receipt.json").read_bytes())
    assert receipt["data_sha256"]==hash_bytes(data)
    assert receipt["pid"]==result["pid"]
    assert receipt["exited"] and receipt["compile_exit_code"]==receipt["run_exit_code"]==0
    source=Path(__file__).with_name("R50FourLineMetal.swift")
    runner=Path(__file__).with_name("run_shard.py")
    assert receipt["source_sha256"]==hash_bytes(source.read_bytes())
    assert receipt["runner_sha256"]==hash_bytes(runner.read_bytes())
    assert result["all_checks_passed"] and result["gpu_command_buffers"]>0
    assert result["gpu_command_seconds_sum"]>0
    for witness in result["witnesses"].values():
        verify_witness(witness)
    for witness in result["direct_cpu_full_walsh_sample_evidence"]:
        coefficients=walsh(table(witness["case_id"]))
        assert len(coefficients)==witness["coefficients"]==128
        assert hash_bytes(struct.pack("<128i",*coefficients))==witness["spectrum_sha256_le_int32"]
    lo,hi=result["case_interval_half_open"]
    assert hi-lo==result["verified_cases"]
    assert result["cpu_gpu_integer_coefficients_verified"]==128*(hi-lo)
    assert result["tested_bands_including_empty"]==16*(hi-lo)
    assert sum(result["walsh_sparsity_histogram"])==hi-lo
    assert result["private_bands_with_two_local_preimages"]<=result["private_bands"]
    assert result["bands_with_complete_coefficient_cancellation"]<=result["colliding_bands"]
    return result,receipt


def main():
    p=argparse.ArgumentParser()
    p.add_argument("local",type=Path)
    p.add_argument("macbook",type=Path)
    p.add_argument("out",type=Path)
    args=p.parse_args()
    a,x=read(args.local)
    b,y=read(args.macbook)
    split=tuple(sorted([a["case_interval_half_open"],b["case_interval_half_open"]]))
    assert split==([0,786432],[786432,1048576])
    assert x["start_epoch"]==y["start_epoch"]
    assert a["bank_truth_sha256_le_int32"]==b["bank_truth_sha256_le_int32"]
    fields=["verified_cases","cpu_gpu_integer_coefficients_verified","tested_bands_including_empty",
            "tested_nonempty_bands","private_bands","private_bands_with_two_local_preimages",
            "colliding_bands","bands_with_complete_coefficient_cancellation",
            "completely_canceled_coefficient_occurrences","band_line_mass_equalities",
            "band_actual_address_mass_equalities","global_line_mass_equalities","direct_cpu_full_walsh_cases"]
    overlap=[(r["actual_start_unix"],r["actual_end_unix"]) for r in (x,y)]
    common=min(t[1] for t in overlap)-max(t[0] for t in overlap)
    out={
        "schema":"a28-r050-verified-two-machine-four-line-sweep-v1",
        "case_intervals":split,"bank_truth_sha256":a["bank_truth_sha256_le_int32"],
        "scheduled_unix":x["start_epoch"],
        "nominal_process_overlap_seconds":max(0,common),
        "nominal_metal_command_overlap_seconds":intersection_time(
            merge(a["gpu_command_utc_unix_intervals"]),merge(b["gpu_command_utc_unix_intervals"])),
        "combined":{field:a[field]+b[field] for field in fields},
        "machine":{result["device"]:{
            "case_interval":result["case_interval_half_open"],"process":receipt["pid"],
            "wall_seconds":receipt["actual_end_unix"]-receipt["actual_start_unix"],
            "gpu_command_buffer_count":result["gpu_command_buffers"],
            "gpu_command_seconds_summed":result["gpu_command_seconds_sum"],
            "python_independent_witness_cases":len(result["direct_cpu_full_walsh_sample_evidence"]),
            "exact_bands":result["witnesses"],"results_sha256":receipt["data_sha256"],
            "process_exited":receipt["exited"],
        } for result,receipt in ((a,x),(b,y))},
    }
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps(out,sort_keys=True))


if __name__=="__main__":
    main()
