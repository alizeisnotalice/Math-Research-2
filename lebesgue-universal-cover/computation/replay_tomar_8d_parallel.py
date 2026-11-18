"""Run every original 8D Arb part with a bounded process pool and replay the exact tree."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("checkout",type=Path)
parser.add_argument("certdir",type=Path)
parser.add_argument("receipts",type=Path)
parser.add_argument("--workers",type=int,default=6)
args=parser.parse_args()
assert args.workers>=1
directory=args.certdir.resolve()
scripts=args.checkout.resolve()/"checker"
args.receipts.mkdir(parents=True,exist_ok=True)
env={**os.environ,"VECLIB_MAXIMUM_THREADS":"1","OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"}
groups=sorted(json.loads(p.read_bytes())["part"] for p in directory.glob("part*.json"))
assert groups==list(range(24))
for name in ["treecheck.py","lbcheck.py"]:
    fixed=subprocess.check_output(["git","-C",str(args.checkout),"show","10af98606c98ca9ffd0d9f402cae1e2efe74b792:checker/"+name])
    assert fixed==(scripts/name).read_bytes()
expected="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
archive=args.checkout/"certificates/cert8d_0836.tar.xz"
archive_sha=hashlib.sha256(archive.read_bytes()).hexdigest()
archive_manifest=subprocess.check_output(["git","-C",str(args.checkout),"show",expected+":certificates/SHA256SUMS"]).decode().splitlines()
assert any(line.split()==[archive_sha,archive.name] for line in archive_manifest)
manifest=subprocess.check_output(["tar","-xOf",str(archive),"cert8d_0836/SHA256SUMS"])
assert manifest==(directory/"SHA256SUMS").read_bytes()
for line in manifest.decode().splitlines():
    sha,name=line.split()
    assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==sha

def run(name,*extra):
    command=[sys.executable,"-I","-B",str(scripts/name),str(directory),*extra]
    start=time.monotonic()
    done=subprocess.run(command,text=True,capture_output=True,check=True,env=env)
    lines=done.stdout.splitlines()
    assert len(lines)==1
    return json.loads(lines[0]),round(time.monotonic()-start,3)

def check_part(p):
    row,seconds=run("lbcheck.py",str(p))
    assert row["part"]==p and row["ok"] and row["fails"]==0
    record={"seconds":seconds,**row}
    final=args.receipts/f"part{p:04}.json"
    tmp=final.with_suffix(".partial")
    tmp.write_text(json.dumps(record,indent=2)+"\n")
    tmp.replace(final)
    print(f"part {p:02} {row['leaves']} leaves {seconds}s",flush=True,file=sys.stderr)
    return record

presence,_=run("lbcheck.py","present")
assert presence[0] is True
start=time.monotonic()
parts_by_work=sorted(groups,key=lambda p:(directory/f"part{p:04}.npy").stat().st_size,reverse=True)
with ThreadPoolExecutor(max_workers=args.workers) as pool:
    records=sorted(pool.map(check_part,parts_by_work),key=lambda x:x["part"])
tree,t_tree=run("treecheck.py")
assert tree["ok"] and tree["header"]["ns"]==[3,5,7] and tree["header"]["tau"]==0.836
assert tree["first"]==0 and tree["last"]==24 and tree["cells"]==65536
assert sum(row["leaves"] for row in records)==tree["leaves"]
import numpy
import flint
print(json.dumps({"frozen_commit":expected,"input_archive_sha256":archive_sha,
    "python":sys.version,"numpy":numpy.__version__,"python_flint":flint.__version__,
    "workers":args.workers,"wall_seconds":round(time.monotonic()-start,3),
    "part_process_elapsed_sum":sum(x["seconds"] for x in records),"tree_seconds":t_tree,
    "parts":records,"tree":tree,"complete":True},indent=2))
