#!/usr/bin/env python3
"""Fresh replay of the degree-nine endpoint sieve and S3 proof."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    if sys.flags.optimize:
        raise SystemExit("Python assertions must be enabled.")
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    if a.workers<1:p.error("workers must be positive")
    out=a.out.resolve()
    if out==ROOT or ROOT in out.parents:
        p.error("choose an output directory outside this package")
    manifest_path=ROOT/"manifest.json"
    manifest=json.loads(manifest_path.read_bytes())
    for name,d in manifest["files_sha256"].items():
        assert digest(ROOT/name)==d, name
    out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    env.update(TMPDIR=str(out),TMP=str(out),TEMP=str(out),PYTHONDONTWRITEBYTECODE="1")
    w=("--workers",str(a.workers))
    jobs=[
        ("envelope","verify_envelope.py","--out",out/"envelope.json"),
        ("singleton-probe","singleton_probe.py","--out",out/"singleton-probe.json"),
        ("witnesses","verify_witnesses.py",out/"singleton-probe.json","--out",out/"witness-verification.json"),
        ("r3-entry","r3_entry.py","--out",out/"r3-entry.json"),
        ("r3-full","r3_full.py","--out",out/"r3-full.json",*w),
        ("r3-trees","r3_trees.py",out/"r3-full.rows.jsonl","--out",out/"r3-trees.json",*w),
        ("r3-independent","verify_r3_independent.py",out/"r3-full.rows.jsonl",out/"r3-trees.rows.jsonl",
         "--out",out/"r3-independent.json",*w),
    ]
    procs=[]
    for label,script,*args in jobs:
        path=ROOT/"source"/script
        cmd=[sys.executable,str(path),*map(str,args)]
        begin=monotonic()
        log,err=out/f"{label}.stdout",out/f"{label}.stderr"
        with log.open("xb") as so,err.open("xb") as se:
            with subprocess.Popen(cmd,env=env,cwd=out,stdout=so,stderr=se) as run:
                code=run.wait()
                pid=run.pid
        if code:
            raise RuntimeError(f"{label} exited {code}. See {err}")
        assert not err.read_bytes(), label
        procs.append(dict(stage=label,pid=pid,exit_code=code,elapsed_seconds=round(monotonic()-begin,6),
                          command=cmd,source_sha256=digest(path)))
    hashes={}
    for stored in sorted((ROOT/"frozen").iterdir()):
        path=out/stored.name
        assert path.read_bytes()==stored.read_bytes(), stored.name
        hashes[stored.name]=digest(path)
    endpoint=json.loads((out/"envelope.json").read_bytes())
    independent=json.loads((out/"r3-independent.json").read_bytes())
    assert len(endpoint["ordinary_supports"])==50 and len(endpoint["reversal_representatives"])==25
    assert independent["first_cores"]==652
    assert independent["counts"]["coupled_mothers"]==1281 and independent["counts"]["root_pass"]==243
    assert independent["independent_tree_deaths"]=={"10":168,"11":48,"12":20,"13":5,"14":2}
    receipt=dict(status="verified",python=sys.version,workers=a.workers,root_survivors=243,
                 supported_orbit_count=25,closed_singleton_orbit=[3,6],source_manifest_sha256=digest(manifest_path),
                 fresh_outputs_sha256=hashes,processes=procs)
    (out/"fresh-receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps(receipt,sort_keys=True))


if __name__=="__main__":
    main()
