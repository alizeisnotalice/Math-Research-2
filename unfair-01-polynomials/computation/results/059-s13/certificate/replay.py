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
FILES = ("four-row-entry.json","four-row-entry.rows.jsonl","independent-entry.json","independent-entry.rows.jsonl")


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
        ("radical", ROOT/"source/s13_exact_four_rows.py", "--out", out/"four-row-entry.json"),
        ("integral-basis", ROOT/"source/verify_s13_integral_basis.py", out/"four-row-entry.json",
         "--out",out/"independent-entry.json")
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
    proof=json.loads((out/"four-row-entry.json").read_bytes())
    witness=json.loads((out/"independent-entry.json").read_bytes())
    assert proof["domains"]==112 and proof["counts"]["mothers"]==0
    assert proof["counts"].get("shared_d",0)==0
    assert (witness["counts"]["shared_unit_d"],witness["counts"]["C4_value_2"])==(866,866)
    receipt=dict(status="verified",python=sys.version,workers=a.workers,domains=proof["domains"],
                 counts=witness["counts"],fresh_outputs_sha256=outputs,
                 source_manifest_sha256=sha(ROOT/"manifest.json"),processes=processes)
    (out/"fresh-receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps(receipt,sort_keys=True))


if __name__=="__main__":
    main()
