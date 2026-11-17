#!/usr/bin/env python3
"""Integrity/coverage verifier for the compact round-055 S={1} archive.

It intentionally does not rerun the prefix or root-screen mathematics.
"""
from __future__ import annotations
import gzip,hashlib,json
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
# Works both from the staging directory and after this file is copied into the
# archive root.
ROOT=HERE if (HERE/"package-manifest.json").is_file() else HERE/"package"
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def lines(rel):
    with gzip.open(ROOT/rel,"rb") as f:
        raw=f.read(); rows=[json.loads(x) for x in raw.splitlines() if x]
    return rows,hashlib.sha256(raw).hexdigest()
def H(row):
    T,N,b,c,d,e,f,g=row; F=[1,0,b,c,d,e,f,g,1]; out=[0]*17
    for i,x in enumerate(F):
        for j,y in enumerate(F):out[i+j]+=x*y
    for i,x in enumerate(F):out[i+1]+=T*x
    out[2]+=N
    return out
def main():
    man=json.loads((ROOT/"package-manifest.json").read_text()); seen=set()
    for x in man["files"]:
        assert x["path"] not in seen;seen.add(x["path"])
        p=ROOT/x["path"];assert p.is_file() and not p.is_symlink() and p.stat().st_size==x["bytes"] and sha(p)==x["sha256"],x["path"]
    assert not any(p.is_symlink() for p in ROOT.rglob("*"))
    cores=json.loads((ROOT/"data/core-list.json").read_text())["cores"]; assert len(cores)==len(set(map(tuple,cores)) )==1140
    c4,h4=lines("data/stage4-core-counts.jsonl.gz"); c57,h57=lines("data/stage57-core-counts.jsonl.gz")
    assert len(c4)==len(c57)==1140 and {tuple(x["core"]) for x in c4}=={tuple(x) for x in cores}=={tuple(x["core"]) for x in c57}
    led4,_=lines("frozen/stage4-original-core-ledger.jsonl.gz"); led57,_=lines("frozen/stage57-original-core-ledger.jsonl.gz")
    def key(x): return (tuple(x["core"]), x.get("records"), x.get("name",x.get("source_name")), x.get("sha256",x.get("source_sha256")), x.get("bytes",x.get("source_bytes")), x["status"])
    assert [key(x) for x in led4]==[key(x) for x in c4] and [key(x) for x in led57]==[key(x) for x in c57], "canonical counts do not bind to original-output ledgers"
    p4=json.loads((ROOT/"data/stage4-parts.json").read_text()); rec4=[]
    for x in p4:
        rr,h=lines(x["file"]); assert h==x["raw_sha256"] and len(rr)==x["cores"]; rec4+=rr
    assert len(rec4)==1140 and sum(len(x["records"]) for x in rec4)==60134 and sum(x["records"] for x in c4)==60134
    idx57=json.loads((ROOT/"data/stage57-records-index.json").read_text());rec57,h=lines(idx57["file"])
    assert h==idx57["raw_sha256"] and len(rec57)==idx57["cores"]==1140 and sum(len(x["records"]) for x in rec57)==idx57["records"]==67 and sum(x["records"] for x in c57)==67
    rindex=json.loads((ROOT/"data/root-mothers-index.json").read_text()); mothers,h=lines(rindex["file"])
    assert h==rindex["raw_sha256"] and len(mothers)==rindex["mothers"]==67
    stage_rows={tuple([x["core"][0],x["core"][1],*r[:6]]) for x in rec57 for r in x["records"]}; assert len(stage_rows)==67
    assert {tuple(x["row"]) for x in mothers}==stage_rows
    codes=Counter()
    for x in mothers:
        assert x["H"]==H(x["row"]) and x["production"]["screen"]["status"]=="rejected"
        codes[x["failure_code"]]+=1
    assert dict(sorted(codes.items()))=={int(k):v for k,v in rindex["codes"].items()}
    dep=json.loads((ROOT/"dependencies/algebraic/full-v1-comparison.json").read_text())
    assert dep["status"]=="passed_all_1140_core_records_and_multiplicities" and dep["row4_records"]==60134 and dep["row7_records"]==67
    dep_cores,_=lines("dependencies/algebraic/full-v1-core-ledger.jsonl.gz")
    assert len(dep_cores)==1140 and sum(x["row4_records"] for x in dep_cores)==60134 and sum(x["row7_records"] for x in dep_cores)==67
    summary=json.loads((ROOT/"archive-summary.json").read_text())
    assert dep["parent_stage4_summary_sha256"]==summary["stage4"]["summary_sha256"]
    assert dep["parent_stage57_summary_sha256"]==summary["stage57"]["summary_sha256"]
    independent=json.loads((ROOT/"frozen/independent-root-summary.json").read_text())
    assert independent["status"]=="passed_all_67" and independent["production_summary_sha256"]==summary["root"]["summary_sha256"]
    assert independent["stage57_sha256"]==summary["stage57"]["summary_sha256"]
    print(json.dumps({"status":"passed","manifest_files":len(man["files"]),"cores":1140,"row4":60134,"row7":67,"root_mothers":67,"code_histogram":dict(sorted(codes.items())),"core_ledger_payloads":{"stage4":h4,"stage57":h57}},sort_keys=True))
if __name__=="__main__":main()
