#!/usr/bin/env python3
"""Bounded R71 common-mass combinations of individually positive r=1/4 cells."""
import json,sys,os,time,hashlib
from pathlib import Path
import joint_outside_lp as j
t=time.monotonic()
left=[[-10,y] for y in [0,24,48,72,96]]
right=[[106,y] for y in [0,24,48,72,96]]
bottom=[[x,-10] for x in [0,24,48,72,96]]
top=[[x,106] for x in [0,24,48,72,96]]
models=[("side-neighbours",left[:2]),("opposite",[left[0],right[0]]),("adjacent-sides",[left[0],bottom[0]]),("side-three",left[:3]),("side-five",left),("opposite-ten",left+right),("adjacent-ten",left+bottom),("perimeter-twenty",left+right+bottom+top)]
results=[]
for name,centres in models:
  r=j.solve(centres,radius=25,export="computation/inputs/r071/"+name if name=="perimeter-twenty" else None);r["name"]=name;results.append(r)
out={"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"module_sha256":hashlib.sha256(Path(j.__file__).read_bytes()).hexdigest(),"python":sys.version,"numpy":j.np.__version__,"scipy":j.scipy.__version__,"pid":os.getpid(),"seconds":time.monotonic()-t,"scope":"8 specified combinations; all r=1/4; each centre individually positive by earlier singleton plus exact D4 symmetries; no search beyond these models", "results":results}
Path("computation/results/round071-joint/shared-batch.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"pid":os.getpid(),"seconds":out["seconds"],"models":[{"name":r["name"],"eta":r.get("eta_float"),"exact_eta":r.get("primal",{}).get("eta_float"),"exact_dual":r.get("dual",{}).get("upper_eta_float"),"centres":len(r["centres_hundredths"])} for r in results]},indent=2))
