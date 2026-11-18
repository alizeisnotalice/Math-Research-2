#!/usr/bin/env python3
"""Complete mod-163 necessary search of the disjoint 65/66 support envelopes."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from itertools import combinations,product
import json
import os
from pathlib import Path
import re
import subprocess
from time import monotonic

ROOT=Path(__file__).resolve().parent


def dump(row):
    return json.dumps(row,sort_keys=True,separators=(",",":")).encode()+b"\n"


def digest(raw):
    return sha256(raw).hexdigest()


def shafile(p):
    return digest(p.read_bytes())


def pack(values):
    return ",".join(map(str,values)) or "-"


def direct_degree(n,support,p=163):
    a=[int(i in support) for i in range(n+1)]
    b=[i*a[i]%p for i in range(1,len(a))]
    while b and b[-1]==0:b.pop()
    while b:
        inverse=pow(b[-1],-1,p)
        while len(a)>=len(b):
            shift=len(a)-len(b)
            factor=a[-1]*inverse%p
            for i,x in enumerate(b):
                a[i+shift]=(a[i+shift]-x*factor)%p
            while a and a[-1]==0:a.pop()
        a,b=b,a
    return len(a)-1


def unpack(cube,gray):
    return sorted(cube["one"]+[i for k,i in enumerate(cube["free"]) if (gray>>k)&1])


def scanner_cmd(binary,n,cube,lo,hi):
    return [str(binary),str(n),pack(cube["one"]),pack(cube["free"]),str(lo),str(hi)]


def check_output(raw,n,cube,lo,hi):
    report=json.loads(raw)
    assert (report["degree"],report["begin"],report["end"],report["count"])==(n,lo,hi,hi-lo)
    ranks=[]
    retained=[]
    for rank,gray,d in report["hits"]:
        assert lo<=rank<hi and gray==rank^(rank>>1)
        ranks.append(rank)
        support=unpack(cube,gray)
        assert d>=7 and d==direct_degree(n,support)
        retained.append(dict(rank=rank,gray=gray,gcd163=d,gcd1279=direct_degree(n,support,1279),support=support))
    assert ranks==sorted(set(ranks))
    return retained


def test(binary,data):
    samples=[]
    for n,steps in ((65,(7,21,37)),(66,(7,21,38))):
        support=sorted(sum(x for x,b in zip(steps,word) if b) for word in product((0,1),repeat=3))
        assert len(support)==len(set(support))
        assert support[-1]==n
        cube=dict(one=support,free=[2,3,9,10,14,15])
        expect=direct_degree(n,support)
        assert expect>=7
        raw=subprocess.check_output(scanner_cmd(binary,n,cube,0,1<<len(cube["free"])))
        got=check_output(raw,n,cube,0,1<<len(cube["free"]))
        expected=[]
        for k in range(1<<len(cube["free"])):
            s=unpack(cube,k^(k>>1))
            d=direct_degree(n,s)
            if d>=7:expected.append((k,d))
        assert [(r["rank"],r["gcd163"]) for r in got]==expected
        assert (0,expect) in expected
        samples.append(dict(n=n,toy_masks=1<<len(cube["free"]),positive_mask_degrees=expected))
    for family in data["degrees"]:
        n=family["degree"]
        checked=0
        for cube in family["disjoint_cubes"]:
            for raw in (0,(1<<len(cube["free"]))-1,0x53627e9&((1<<len(cube["free"]))-1)):
                support=unpack(cube,raw)
                expected=direct_degree(n,support)
                value=subprocess.check_output([str(binary),"--inspect",str(n),pack(support)])
                assert expected==int(value)
                checked+=1
        samples.append(dict(n=n,actual_sample_degrees_verified=checked))
    return samples


def shard(binary,destination,n,k,cube,lo,hi):
    full=destination/f"{n}-{k:03d}-{lo:08d}.json"
    assert not full.exists()
    cmd=["/usr/bin/time","-l",*scanner_cmd(binary,n,cube,lo,hi)]
    t=monotonic()
    with subprocess.Popen(cmd,cwd=destination,stdout=subprocess.PIPE,stderr=subprocess.PIPE) as proc:
        stdout,stderr=proc.communicate(timeout=180)
        pid,code=proc.pid,proc.returncode
    assert code==0, (cmd,stderr)
    hit=check_output(stdout,n,cube,lo,hi)
    m=re.search(rb"(\d+)\s+maximum resident set size",stderr)
    assert m
    result=dict(degree=n,cube=k,begin=lo,end=hi,examined=hi-lo,hits=hit,
                pid=pid,exit_code=code,peak_rss_bytes=int(m[1]),elapsed_seconds=round(monotonic()-t,6),
                stdout_sha256=digest(stdout),command=cmd)
    with full.with_suffix(".partial").open("xb") as f:
        f.write(dump(result))
        f.flush()
        os.fsync(f.fileno())
    full.with_suffix(".partial").replace(full)
    return full,result


def verify_family(group):
    n=group["degree"]
    cs=group["disjoint_cubes"]
    for r in cs:
        a,b,c=(set(r[k]) for k in ("one","zero","free"))
        assert not (a&b or b&c or c&a)
        assert a|b|c==set(range(n+1))
        assert {0,n}<=a and r["count"]==1<<len(c)
    for a,b in combinations(cs,2):
        assert set(a["one"])&set(b["zero"]) or set(b["one"])&set(a["zero"])
    assert sum(r["count"] for r in cs)==group["distinct_products_to_check"]
    assert all(len(r["free"])<=26 for r in cs)


def main():
    if not __debug__:raise SystemExit("assertions required")
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--binary",type=Path,default=ROOT/"sieve")
    ap.add_argument("--workers",type=int,default=3)
    ap.add_argument("--block",type=int,default=1<<19)
    ap.add_argument("--smoke",action="store_true")
    a=ap.parse_args()
    if a.workers<1 or a.block<1:ap.error("workers and block must be positive")
    out=a.out.resolve()
    binary=a.binary.resolve()
    data=json.loads(a.input.read_bytes())
    for family in data["degrees"]:verify_family(family)
    out.mkdir(exist_ok=False)
    (out/"smoke.json").write_bytes(dump(test(binary,data)))
    if a.smoke:
        print((out/"smoke.json").read_text())
        return
    manifest=dict(scope="all distinct products in support envelopes, mod163 gcd>=7 only",
                  support_sha256=shafile(a.input),binary_sha256=shafile(binary),
                  source_sha256=shafile(ROOT/"sieve.cpp"),workers=a.workers,block=a.block)
    (out/"manifest.json").write_bytes(dump(manifest))
    sh=out/"shards"
    sh.mkdir()
    jobs=[(binary,sh,fam["degree"],i,r,lo,min(lo+a.block,r["count"]))
          for fam in data["degrees"] for i,r in enumerate(fam["disjoint_cubes"])
          for lo in range(0,r["count"],a.block)]
    results=[]
    start=monotonic()
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(shard,*job) for job in jobs]
        for i,future in enumerate(futures,1):
            results.append(future.result())
            if i%16==0 or i==len(futures):
                print(json.dumps(dict(ordered_complete=i,total_shards=len(futures),
                                      written_shards=len(list(sh.glob("*.json"))))),flush=True)
    counts=Counter()
    candidates=[]
    index=[]
    pids=[]
    for (path,r),job in zip(results,jobs,strict=True):
        _,_,n,k,_,lo,hi=job
        assert (n,k,lo,hi)==(r["degree"],r["cube"],r["begin"],r["end"])
        counts[n]+=r["examined"]
        candidates.extend(dict(degree=n,cube=k,**row) for row in r["hits"])
        index.append(dict(file=path.name,sha256=shafile(path),stdout_sha256=r["stdout_sha256"],
                          examined=r["examined"],mod163_hits=len(r["hits"])))
        pids.append(r["pid"])
    for pid in pids:
        try:
            os.kill(pid,0)
        except ProcessLookupError:
            pass
        else:
            raise RuntimeError(f"PID {pid} still present")
    assert counts==Counter({fam["degree"]:fam["distinct_products_to_check"] for fam in data["degrees"]})
    (out/"shards-index.jsonl").write_bytes(b"".join(map(dump,index)))
    (out/"hits.jsonl").write_bytes(b"".join(map(dump,candidates)))
    complete=dict(status="complete",candidate_counts=counts,mod163_hits=Counter(r["degree"] for r in candidates),
                  both_prime_hits=Counter(r["degree"] for r in candidates if r["gcd1279"]>=7),
                  worker_processes_exited=len(pids),max_rss_per_child=max(r["peak_rss_bytes"] for _,r in results),
                  elapsed_sieve_seconds=round(monotonic()-start,6),
                  index_sha256=shafile(out/"shards-index.jsonl"),hits_sha256=shafile(out/"hits.jsonl"),
                  manifest_sha256=shafile(out/"manifest.json"))
    (out/"complete.json").write_bytes(dump(complete))
    print(json.dumps(complete,sort_keys=True))


if __name__=="__main__":main()
