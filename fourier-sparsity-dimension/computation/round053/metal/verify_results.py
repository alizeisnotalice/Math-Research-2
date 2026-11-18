"""Independently verify R53 provenance, small-code counts, and real orbit witnesses."""
import argparse
import collections
import hashlib
import itertools
import json
import math
from pathlib import Path

SIZE=1<<14


def fwht(values):
    for step in (1<<k for k in range(14)):
        for start in range(0,len(values),2*step):
            for x in range(start,start+step):
                a,b=values[x],values[x+step]
                values[x],values[x+step]=a+b,a-b
    return values


def echelon(values):
    result={}
    for value in values:
        while value:
            top=value.bit_length()-1
            if top in result:
                value^=result[top]
            else:
                result[top]=value
                break
    return result


def affine_rank(points):
    anchor=points[0]
    return len(echelon(v^anchor for v in points))


def truth(row):
    seed=[0]*SIZE
    for s,c in zip(row["support"],row["signs"]):
        seed[s>>1]=c
    values=fwht(seed)
    assert all(abs(w)==8 for w in values)
    return [int(a<0) for a in values]


def signed_words(row):
    words=[0]
    for j in range(15):
        col=sum(((row["support"][x]>>j)&1)<<x for x in range(64))
        words.extend([word^col for word in words])
    assert len(set(words))==32768
    min_words=[i for i,word in enumerate(words) if word.bit_count()==8]
    assert min(w.bit_count() for w in words[1:])==8
    return words,min_words


def phase_activity(table, direction):
    a,c=direction>>1,direction&1
    result=bytearray(SIZE//8)
    active=[]
    for x in range(SIZE):
        if table[x]^table[x^a]^c:
            result[x>>3]|=1<<(x&7)
            active.append(x)
    return int.from_bytes(result,"little"),active


def span(basis):
    out=[0]
    for d in basis:
        out += [x^d for x in out]
    return out


def check_witness(row,words,table,witness):
    basis=witness["minimum_phase_direction_basis"]
    cs=span(basis)
    assert all(words[d].bit_count()==8 for d in basis)
    assert cs[1:]==witness["nonzero_directions"]
    assert [words[d].bit_count() for d in cs[1:]]==witness["weights"]
    union=0
    for d in cs:
        union|=words[d]
    assert union.bit_count()==witness["spectral_support_union"]
    pivots=echelon(d>>1 for d in basis)
    skip=sum(1<<i for i in pivots)
    assert len(pivots)==len(basis)
    counter=collections.Counter()
    for x in range(SIZE):
        if x&skip: continue
        orbit=sum((table[x^(d>>1)]^(d&1))<<i for i,d in enumerate(cs))
        counter[str(orbit)]+=1
    assert dict(counter)==witness["twisted_truth_pattern_histogram"]
    assert sum(counter.values())==witness["physical_orbit_count"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(directory, source, driver, seeds):
    data=directory/"result.json"
    result=json.loads(data.read_bytes())
    receipt=json.loads((directory/"receipt.json").read_bytes())
    assert receipt["data_sha256"]==digest(data)
    assert receipt["run_exit_code"]==receipt["compile_exit_code"]==0 and receipt["exited"]
    assert receipt["source_sha256"]==digest(source)
    assert receipt["runner_sha256"]==digest(driver)
    assert receipt["input_sha256"]==result["input_sha256"]==digest(seeds)
    assert result["all_checks_passed"] and result["gpu_orbit_command_seconds"]>0
    return result,receipt


def main():
    parser=argparse.ArgumentParser()
    for name in ("local","macbook","out"):
        parser.add_argument(name,type=Path)
    opts=parser.parse_args()
    root=Path(__file__).resolve().parent
    seeds=root/"r11_minimal_input.json"
    read_args=(root/"R53FlatOrbitMetal.swift",root/"run_shard.py",seeds)
    a,ar=read(opts.local,*read_args)
    b,br=read(opts.macbook,*read_args)
    source=json.loads(seeds.read_bytes())
    assert source["source_sha256"]==a["archive_input_sha256"]==b["archive_input_sha256"]
    assert a["row_interval_half_open"]==[15,20] and b["row_interval_half_open"]==[0,15]
    assert ar["start_epoch"]==br["start_epoch"]
    input_by_id={row["index"]:row for row in source["rows"]}
    output_by_id={row["row_index"]:row for host in (a,b) for row in host["rows"]}
    assert input_by_id.keys()==output_by_id.keys()==set(range(20))
    pair_unions=collections.Counter()
    triple_unions=collections.Counter()
    pair_intersections=collections.Counter()
    triple_intersections=collections.Counter()
    hist_cache={}
    phase_counts=collections.Counter()
    for i,row in input_by_id.items():
        words,min_words=signed_words(row)
        got=output_by_id[i]
        assert got["minimum_phase_directions"]==min_words and len(min_words)==24
        assert got["period_dimension"]==row["period_dimension"]
        assert dict(collections.Counter(str(w.bit_count()) for w in words))==got["full_direction_weight_histogram"]
        table=truth(row)
        hist_cache[i]=(row,words,table)
        min_masks={}
        for direction in min_words:
            active,points=phase_activity(table,direction)
            assert active.bit_count()==len(points)==2048 and affine_rank(points)==11
            frequencies=[s>>1 for j,s in enumerate(row["support"]) if (words[direction]>>j)&1]
            assert len(frequencies)==8 and affine_rank(frequencies)==3
            phase_counts[direction&1]+=1
            min_masks[direction]=active
        rank2=set()
        rank3=set()
        for d,e in itertools.combinations(min_words,2):
            rank2.add(tuple(sorted((d,e,d^e))))
            pair_intersections[(min_masks[d]&min_masks[e]).bit_count()]+=1
        for d,e,f in itertools.combinations(min_words,3):
            if d^e==f: continue
            rank3.add(tuple(sorted((d,e,f,d^e,d^f,e^f,d^e^f))))
            triple_intersections[(min_masks[d]&min_masks[e]&min_masks[f]).bit_count()]+=1
        assert [len(min_words),len(rank2),len(rank3)]==got["row_codes_dim1_dim2_dim3"]
        for counter,spaces in ((pair_unions,rank2),(triple_unions,rank3)):
            for code in spaces:
                union=0
                for d in code:
                    union|=words[d]
                counter[union.bit_count()]+=1
    for rank,expected in ((2,pair_unions),(3,triple_unions)):
        for weight,count in expected.items():
            key=f"dim{rank}_union{weight}"
            assert a["orbit_statistics_by_dimension_and_spectral_union"][key]["subcodes"] +                    b["orbit_statistics_by_dimension_and_spectral_union"][key]["subcodes"] == count
    for host in (a,b):
        for witness in host["orbit_exact_witnesses"].values():
            check_witness(*hist_cache[witness["row"]],witness)
    for phase,n in phase_counts.items():
        assert a["minimum_directions_by_phase_bit"][phase]+b["minimum_directions_by_phase_bit"][phase]==n
    totals={}
    for s in ("complete_subcode_counts_dim1_dim2_dim3","metal_orbits_dim1_dim2_dim3",
              "metal_twisted_truth_values_dim1_dim2_dim3","direct_cpu_orbit_histograms_crosschecked_dim1_dim2_dim3"):
        totals[s]=[x+y for x,y in zip(a[s],b[s])]
    for s in ("whole_boolean_truth_cpu_metal_values","enumerated_phase_translation_words_including_zero",
              "minimum_spectral_flats_and_physical_flats_checked","derivative_density_equalities_all_nonzero_code_directions"):
        totals[s]=a[s]+b[s]
    result={
        "schema":"a28-r053-independently-verified-phase-orbits-v1",
        "archive_input_sha256":source["source_sha256"],"totals":totals,
        "minimum_directions_by_phase_bit":dict(phase_counts),
        "rank2_spectral_union_distribution":dict(sorted(pair_unions.items())),
        "rank3_spectral_union_distribution":dict(sorted(triple_unions.items())),
        "all_minimum_direction_pairs_by_physical_coactivity_point_count":dict(sorted(pair_intersections.items())),
        "all_independent_minimum_direction_triples_by_physical_coactivity_point_count":dict(sorted(triple_intersections.items())),
        "python_cpu_metal_witness_histograms_checked":sum(len(host["orbit_exact_witnesses"]) for host in (a,b)),
        "nominal_process_overlap_seconds":
            min(ar["actual_end_unix"],br["actual_end_unix"])-max(ar["actual_start_unix"],br["actual_start_unix"]),
        "devices":{host["device"]:{"rows":host["row_interval_half_open"],"pid":r["pid"],
                   "fwht_gpu_seconds":host["gpu_fwht_command_seconds"],"orbital_gpu_seconds":host["gpu_orbit_command_seconds"],
                   "process_exited":r["exited"],"raw_output_sha256":r["data_sha256"]} for host,r in ((a,ar),(b,br))},
    }
    opts.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
