"""Recompute the independent five-body 0.8372 universal lower bound from the original trees."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
for k in ("cert","slice_index","work","cache"):p.add_argument(k,type=Path)
p.add_argument("--workers",type=int,default=2)
a=p.parse_args();a.work.mkdir(exist_ok=False);w=a.work.resolve();here=Path(__file__).resolve().parent
def run(script,*args,log):
    file=w/log
    with file.with_suffix(".partial").open("w") as out:
        subprocess.run([sys.executable,"-I","-B",here/script,*map(str,args)],stdout=out,check=True)
    file.with_suffix(".partial").replace(file)
low=w/"lower"
run("replay_tomar_grouped_full.py",a.cert,a.slice_index,low,a.cache,"--workers",a.workers,log="lower-replay.json")
parent=low/"audit.json";native=low/"native";table=low/"witnesses32.npz";primary=w/"primary";extra=w/"extra"
record=w/"primary-ledger.json";selection=w/"selection.json";archive=w/"tomar-grouped-08372.zip"
run("lift_tomar_grouped_scores.py",parent,native,a.cert,table,primary,a.cache,"--target","0.8372","--calls",1048575,"--depth",60,log="primary-run.log")
run("audit_tomar_score_lift.py",parent,native,a.cert,primary,w,a.cache,"--target","0.8372","--partial",log=record.name)
run("select_tomar_score_exceptions.py",parent,record,primary,selection,log="selection.log")
run("check_tomar_score_exceptions.py",selection,a.cert,table,a.slice_index,here/"tomar_fine_angle_fan.cpp",extra,a.cache,log="exceptions.log")
run("seal_tomar_final_lift.py",record,primary,selection,extra,archive,log="sealed.json")
run("audit_tomar_final_lift.py",a.cert,table,a.slice_index,parent,low/"public/tomar-grouped-08371-index.json",archive,w,a.cache,
    "--replay-exceptions",log="final.json")
print((w/"final.json").read_text())
