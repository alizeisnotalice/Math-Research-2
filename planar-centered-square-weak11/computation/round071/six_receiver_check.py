#!/usr/bin/env python3
"""R71 compression to dual-active receivers and six prescribed deletions."""
import json,sys,os,time,hashlib
from pathlib import Path
import joint_outside_lp as j
t=time.monotonic()
centres=[[-10,24],[-10,48],[-10,96],[0,106],[48,106],[72,106]]
models=[("six-active",centres)]+[("omit-"+str(i),centres[:i]+centres[i+1:]) for i in range(6)]
results=[]
for name,ps in models:
  r=j.solve(ps,radius=25,export="computation/inputs/r071/six-active" if name=="six-active" else None);r["name"]=name;results.append(r)
out={"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"module_sha256":hashlib.sha256(Path(j.__file__).read_bytes()).hexdigest(),"python":sys.version,"numpy":j.np.__version__,"scipy":j.scipy.__version__,"pid":os.getpid(),"seconds":time.monotonic()-t,"scope":"one dual-selected six-centre model and all six one-centre omissions, r=1/4 fixed", "results":results}
Path("computation/results/round071-joint/six-receiver-check.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"pid":os.getpid(),"seconds":out["seconds"],"models":[{"name":r["name"],"eta":r.get("eta_float"),"exact_eta":r.get("primal",{}).get("eta_float"),"exact_dual":r.get("dual",{}).get("upper_eta_float")} for r in results]},indent=2))
