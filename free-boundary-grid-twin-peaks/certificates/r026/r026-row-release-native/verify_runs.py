#!/usr/bin/env python3
import json
from pathlib import Path

known={3:("156","6"),4:("4213440","5600"),5:("153450311313600","3971406600")}
out=[]
for n,(ht,hr) in known.items():
    a=json.loads(Path(f"n{n}-t1.json").read_text())
    b=json.loads(Path(f"n{n}-t2.json").read_text())
    fields=("n","M","algorithm","states_visited","H_two_rows","H_toprow","limit30")
    assert {k:a[k] for k in fields}=={k:b[k] for k in fields}
    assert (a["H_toprow"],a["H_two_rows"])==(ht,hr)
    out.append({"n":n,"T1_T2_match":True,"H_toprow":ht,"H_two_rows":hr,"states_visited":a["states_visited"]})
print(json.dumps({"all_match":True,"runs":out},indent=2))
