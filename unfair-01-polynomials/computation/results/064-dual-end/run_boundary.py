"""Run each admissible degree-64 parameter once, saving completed output."""
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import time

import sympy as sp

root=Path(__file__).resolve().parent
binary=root/"search_j2"
source=root/"search_j2_degree64.cpp"
out=root/"out"
out.mkdir()
x=sp.Symbol("x")
q_values=[]
for q in range(163):
    j=sp.Poly.from_list([1,6,121,106,68,115,117,q],x,modulus=163)
    if j.nth(0) and j.gcd(j.diff()).degree()==0:
        q_values.append(q)
assert len(q_values)==162
manifest={"degree":64,"prime":163,"q_values":q_values,"left_per_q":1<<26,"right_per_q":1<<26,
          "source_sha256":sha256(source.read_bytes()).hexdigest(),
          "binary_sha256":sha256(binary.read_bytes()).hexdigest(),
          "scope":"Squarefree in Q, boundary gcd degree 7 at 163. Stage tests J^2 only."}
(root/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
for q in q_values:
    start=time.monotonic()
    # Each child owns its table and releases its memory on exit.
    result=subprocess.run([str(binary),"64",str(q)],check=True,capture_output=True,text=True,timeout=14400)
    row=json.loads(result.stdout)
    assert row["degree"]==64 and row["q"]==q and row["mode"]=="boundary"
    assert row["left_patterns"]==row["right_patterns"]==1<<26
    assert len(row["matches"])==len(set(row["matches"]))
    assert all(0<=m<1<<52 for m in row["matches"])
    row["seconds"]=time.monotonic()-start
    path=out/f"q{q:03d}.json"
    with path.with_suffix(".partial").open("x") as f:
        json.dump(row,f,indent=2)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    path.with_suffix(".partial").replace(path)
    print(json.dumps({"q":q,"matches":len(row["matches"]),"seconds":row["seconds"]}),flush=True)
(root/"complete.json").write_text(json.dumps({"status":"full_boundary_stage","parameters":len(q_values)})+"\n")
