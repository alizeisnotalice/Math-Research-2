#!/usr/bin/env python3
"""Recompute both complete covers from the source in this portable package."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
FILES = ("entry.json","entry.rows.jsonl","trees.json","trees.rows.jsonl",
         "independent-entry.json","independent-trees.json")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if sys.flags.optimize:
        raise SystemExit("Python assertions must be enabled.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers",type=int,default=8)
    a = ap.parse_args()
    if a.workers<1:
        ap.error("workers must be positive")
    out = a.out.resolve()
    if out==ROOT or ROOT in out.parents:
        ap.error("choose an output directory outside the package")
    baseline = json.loads((ROOT/"manifest.json").read_bytes())
    for name,digest in baseline["files_sha256"].items():
        assert sha(ROOT/name)==digest, name
    out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    env.update(TMPDIR=str(out),TMP=str(out),TEMP=str(out),PYTHONDONTWRITEBYTECODE="1")
    options=["--workers",str(a.workers)]
    commands = [
        ("entry",ROOT/"source/s2345_entry.py","--out",out/"entry.json"),
        ("trees",ROOT/"source/s2345_trees.py",out/"entry.rows.jsonl","--out",out/"trees.json","--depth","128"),
        ("independent-entry",ROOT/"source/verify_s2345_entry.py",out/"entry.rows.jsonl","--out",
         out/"independent-entry.json"),
        ("independent-trees",ROOT/"source/verify_s2345_trees.py",out/"trees.rows.jsonl","--out",
         out/"independent-trees.json"),
    ]
    processes=[]
    for label,source,*args in commands:
        command=[sys.executable,str(source),*map(str,args),*options]
        start=time.monotonic()
        with (out/f"{label}.stdout").open("xb") as stdout, (out/f"{label}.stderr").open("xb") as stderr:
            with subprocess.Popen(command,env=env,cwd=out,stdout=stdout,stderr=stderr) as proc:
                pid=proc.pid
                code=proc.wait()
        result=dict(stage=label,pid=pid,exit_code=code,elapsed_seconds=round(time.monotonic()-start,6),
                    source_sha256=sha(source),command=command)
        processes.append(result)
        if code:
            raise RuntimeError(f"{label} exited {code}")
        assert not (out/f"{label}.stderr").read_bytes(), label
    outputs={}
    for name in FILES:
        p=out/name
        assert p.read_bytes()==(ROOT/"frozen"/name).read_bytes(), name
        outputs[name]=sha(p)
    entry=json.loads((out/"entry.json").read_bytes())
    trees=json.loads((out/"trees.json").read_bytes())
    independently=json.loads((out/"independent-trees.json").read_bytes())
    assert entry["domains"]==88 and entry["counts"]["a_cores"]==11080
    assert entry["counts"]["mothers"]==421 and entry["counts"]["root_pass"]==101
    assert trees["counts"]=={"dead":101,"forward_dead":101,"reverse_dead":101}
    assert independently["rows"]==101 and independently["status"]=="verified"
    assert max(int(k) for k in independently["forward"])==92
    assert max(int(k) for k in independently["reverse"])==28
    receipt=dict(status="verified",python=sys.version,workers=a.workers,domains=entry["domains"],
                 counts=entry["counts"],trees=trees["counts"],fresh_outputs_sha256=outputs,
                 source_manifest_sha256=sha(ROOT/"manifest.json"),processes=processes)
    (out/"fresh-receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps(receipt,sort_keys=True))


if __name__=="__main__":
    main()
