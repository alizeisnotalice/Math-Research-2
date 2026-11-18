#!/usr/bin/env python3
"""Authenticate full coverage and deterministic receipts without redoing the sieve."""

import argparse
from collections import Counter
import json
from pathlib import Path

from generate_support import one_degree
from run_sieve import dump,pack,shafile,digest,verify_family,check_output


def main():
    if not __debug__:raise SystemExit("assertions required")
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--run",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--export-stable",type=Path)
    p.add_argument("--compare-stable",type=Path)
    a=p.parse_args()
    assert not a.out.exists()
    original=json.loads(a.input.read_bytes())
    regenerated=json.loads(dump(dict(status="verified",degrees=[one_degree(n) for n in (65,66)])))
    assert original==regenerated
    manifest=json.loads((a.run/"manifest.json").read_bytes())
    summary=json.loads((a.run/"complete.json").read_bytes())
    assert manifest["support_sha256"]==shafile(a.input)
    assert manifest["source_sha256"]==shafile(Path(__file__).with_name("sieve.cpp"))
    assert summary["manifest_sha256"]==shafile(a.run/"manifest.json")
    assert summary["index_sha256"]==shafile(a.run/"shards-index.jsonl")
    assert summary["hits_sha256"]==shafile(a.run/"hits.jsonl")
    block=manifest["block"]
    paths=a.run/"shards"
    jobs=[]
    for family in original["degrees"]:
        verify_family(family)
        for index,cube in enumerate(family["disjoint_cubes"]):
            for lo in range(0,cube["count"],block):
                jobs.append((family["degree"],index,cube,lo,min(lo+block,cube["count"])))
    rows=list(map(json.loads,(a.run/"shards-index.jsonl").read_bytes().splitlines()))
    stable=[]
    seen=Counter()
    all_hits=[]
    assert len(rows)==len(jobs)==len(list(paths.glob("*")))==summary["worker_processes_exited"]
    for row,(degree,k,cube,lo,hi) in zip(rows,jobs,strict=True):
        name=f"{degree}-{k:03d}-{lo:08d}.json"
        raw=paths/name
        assert row["file"]==name and row["sha256"]==shafile(raw)
        shard=json.loads(raw.read_bytes())
        assert shard["exit_code"]==0
        assert shard["command"][-5:]==[str(degree),pack(cube["one"]),pack(cube["free"]),str(lo),str(hi)]
        rec=dict(degree=degree,begin=lo,end=hi,
                 hits=[[d["rank"],d["gray"],d["gcd163"]] for d in shard["hits"]],count=hi-lo)
        stdout=json.dumps(rec,separators=(",",":")).encode()+b"\n"
        assert check_output(stdout,degree,cube,lo,hi)==shard["hits"]
        assert digest(stdout)==shard["stdout_sha256"]==row["stdout_sha256"]
        assert shard["examined"]==row["examined"]==hi-lo
        assert len(shard["hits"])==row["mod163_hits"]
        seen[str(degree)]+=hi-lo
        all_hits.extend(dict(degree=degree,cube=k,**x) for x in shard["hits"])
        stable.append({k:v for k,v in row.items() if k!="sha256"})
    assert seen==summary["candidate_counts"]
    assert all_hits==list(map(json.loads,(a.run/"hits.jsonl").read_bytes().splitlines()))
    assert not all_hits and not summary["mod163_hits"]
    canonical=b"".join(map(dump,stable))
    if a.export_stable:
        with a.export_stable.open("xb") as f:f.write(canonical)
    if a.compare_stable:
        assert canonical==a.compare_stable.read_bytes()
    result=dict(status="verified",candidates=dict(seen),hits=len(all_hits),shards=len(rows),
                stable_sha256=digest(canonical),index_sha256=summary["index_sha256"])
    with a.out.open("xb") as f:f.write(dump(result))
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
