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
FILES = ("front.json","front.rows.jsonl","tail.json","tail.rows.jsonl","trees.json","trees.rows.jsonl",
         "independent-front.json","independent-tail.json","independent-trees.json")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if sys.flags.optimize:
        raise SystemExit("Python assertions must be enabled.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers",type=int,default=4)
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
        ("front",ROOT/"source/s23_entry.py","--out",out/"front.json"),
        ("tail",ROOT/"source/s23_tail.py",out/"front.rows.jsonl","--out",out/"tail.json"),
        ("trees",ROOT/"source/s23_trees.py",out/"tail.rows.jsonl","--out",out/"trees.json"),
        ("independent-front",ROOT/"source/verify_s23_front.py",out/"front.rows.jsonl","--out",
         out/"independent-front.json"),
        ("independent-tail",ROOT/"source/verify_s23_tail.py",out/"front.rows.jsonl",out/"tail.rows.jsonl",
         "--out",out/"independent-tail.json"),
        ("independent-trees",ROOT/"source/verify_s23_trees.py",out/"trees.rows.jsonl","--out",
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
    front=json.loads((out/"front.json").read_bytes())
    tail=json.loads((out/"tail.json").read_bytes())
    closed=json.loads((out/"independent-trees.json").read_bytes())
    assert front["domains"]==88 and front["a_cores"]==11760 and front["counts"]["row5"]==3980
    assert tail["mothers"]==5203 and tail["codes"]["255"]==849
    assert closed["counts"]=={"dead_forward_or_reverse":848,"dead_after_factor":1}
    assert closed["exceptional"]["tree"]["depth"]==67
    receipt=dict(status="verified",python=sys.version,workers=a.workers,domains=front["domains"],
                 counts=closed["counts"],fresh_outputs_sha256=outputs,
                 source_manifest_sha256=sha(ROOT/"manifest.json"),processes=processes)
    (out/"fresh-receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps(receipt,sort_keys=True))


if __name__=="__main__":
    main()
