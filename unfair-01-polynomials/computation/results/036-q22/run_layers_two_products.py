#!/usr/bin/env python3
"""Exact q-layer support ledger plus parent seven-rule discovery runs.

The parent program is deliberately external and is used only to discover
which necessary-rule traces remain unknown.  This program owns the DFS input
coverage, canonical hashes, checkpoints, and stop-after-unknown policy.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from math import gcd
from pathlib import Path

PARENT = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round035-parent/propagate_units_two_products.py")


def digest_entries(entries):
    payload = "".join(json.dumps(x, sort_keys=True, separators=(",", ":")) + "\n" for x in entries)
    return hashlib.sha256(payload.encode()).hexdigest()


def dfs_T(p, q):
    selected = {0}; diffs = {p, q-p}
    def walk(x):
        if x == q:
            yield tuple(sorted(selected)); return
        if all(x-d not in selected for d in diffs):
            selected.add(x); yield from walk(x+1); selected.remove(x)
        yield from walk(x+1)
    yield from walk(1)


def G_of(T, q): return tuple(sorted(set(T) | {x+q for x in T}))


def no_singleton(S, p):
    counts = [0]*p
    for x in S: counts[x % p] += 1
    return all(v != 1 for v in counts)


def enumerate_layer(q):
    rows=[]; canonical=[]; total=0
    for p in range(2, q):
        if gcd(p,q) != 1 or q < 2*p: continue
        survivors=[]
        for T in dfs_T(p,q):
            total += 1; S=G_of(T,q)
            if no_singleton(S,p):
                item={"T":list(T),"G_support":list(S)}
                survivors.append(item); canonical.append({"p":p,"q":q,**item})
        rows.append({"p":p,"q":q,"no_mod_p_singleton_survivors":survivors})
    return rows,canonical,total


def dump(path, data): path.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")


def run_one(root, q):
    qdir=root/f"q{q:02d}"; assert not qdir.exists(), qdir
    qdir.mkdir(); started=time.time()
    rows,canonical,total=enumerate_layer(q)
    input_data={"schema":"round036-q-layer-input-v1","q":q,"rows":rows}
    input_path=qdir/"input.json"; dump(input_path,input_data); raw=input_path.read_bytes()
    checkpoint={"status":"enumeration_complete","q":q,"total_T":total,"survivors":len(canonical),"canonical_stream_sha256":digest_entries(canonical),"input_sha256":hashlib.sha256(raw).hexdigest(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"python":sys.version,"started_unix":started}
    dump(qdir/"checkpoint.json",checkpoint)
    out=qdir/"parent-two-products.json.gz"
    cmd=[sys.executable,str(PARENT),"--input",str(input_path.resolve()),"--output",str(out.resolve())]
    proc=subprocess.run(cmd,cwd=qdir,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    (qdir/"parent.stdout.txt").write_text(proc.stdout); (qdir/"parent.stderr.txt").write_text(proc.stderr)
    checkpoint.update({"status":"propagation_complete" if proc.returncode==0 else "propagation_failed","parent_command":cmd,"parent_returncode":proc.returncode,"parent_output_sha256":hashlib.sha256(out.read_bytes()).hexdigest() if out.exists() else None,"elapsed_seconds":time.time()-started})
    dump(qdir/"checkpoint.json",checkpoint)
    if proc.returncode: raise RuntimeError(f"parent failed for q={q}")
    parent_data=json.loads(gzip.decompress(out.read_bytes()))
    lookup={(r["p"],tuple(item["G_support"])):item["T"] for r in rows for item in r["no_mod_p_singleton_survivors"]}
    unknown=next((r for r in parent_data["rows"] if r["unknown_indices"]),None)
    first=None if unknown is None else {"p":unknown["p"],"q":unknown["q"],"T":lookup[(unknown["p"],tuple(unknown["G_support"]))],"G_support":unknown["G_support"],"unit_indices":unknown["unit_indices"],"unknown_indices":unknown["unknown_indices"]}
    result={"q":q,"total_T":total,"survivors":len(canonical),"canonical_stream_sha256":digest_entries(canonical),"parent_all_unit":parent_data["all_unit"],"parent_unknown":parent_data["unknown"],"first_unknown":first,"checkpoint":str((qdir/"checkpoint.json").resolve())}
    dump(qdir/"layer-summary.json",result)
    return result


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,required=True); ap.add_argument("--q-start",type=int,default=22); ap.add_argument("--q-end",type=int,default=26); args=ap.parse_args()
    assert Path("/Volumes/data").is_mount() and PARENT.is_file() and 3<=args.q_start<=args.q_end<=26
    summaries=[]
    for q in range(args.q_start,args.q_end+1):
        summary=run_one(args.root,q); summaries.append(summary)
        if summary["parent_unknown"]:
            break
    manifest={"schema":"round036-two-products-manifest-v1","q_start":args.q_start,"q_end":args.q_end,"completed_qs":[x["q"] for x in summaries],"stopped_on_unknown":bool(summaries and summaries[-1]["parent_unknown"]),"layers":summaries,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dump(args.root/"manifest.json",manifest)
    print(json.dumps({"completed_qs":manifest["completed_qs"],"stopped_on_unknown":manifest["stopped_on_unknown"],"last":summaries[-1] if summaries else None},sort_keys=True))


if __name__=="__main__": main()
