"""Replay every accepted exception tree with a bounded set of independent workers."""
import argparse,hashlib,json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,sys

p=argparse.ArgumentParser(description=__doc__)
for k in ("checkout","cert","table","high_dir","archive","scratch","cache"):p.add_argument(k,type=Path)
p.add_argument("--workers",type=int,default=2)
a=p.parse_args();assert 1<=a.workers<=8
raw=json.loads((a.high_dir/"five-global-localization.json").read_bytes())["retained"]
groups=[[] for _ in range(a.workers)];loads=[0]*a.workers
for part,ids in sorted(raw.items(),key=lambda kv:-len(kv[1])):
    j=min(range(a.workers),key=lambda n:loads[n]);groups[j].append(int(part));loads[j]+=len(ids)
args=[a.checkout,a.cert,a.table,a.high_dir,a.archive,a.scratch,a.cache]
def task(ids):
    r=json.loads(subprocess.check_output([sys.executable,"-I","-B",Path(__file__).with_name("audit_tomar_five_083729.py"),*args,
        "--replay-parts",",".join(map(str,ids))]))
    assert r["complete"] and r["replayed_exception_root_parts"]==sorted(ids)
    return r
with ThreadPoolExecutor(max_workers=a.workers)as pool:results=list(pool.map(task,groups))
assert sum(r["successful_original_trees_replayed"]for r in results)==501
assert len({r["archive_sha256"]for r in results})==1
print(json.dumps({"target":"0.83729","complete":True,"all_exception_arithmetic_replayed":True,
    "high_target_parent_replayed_in_parent_receipt":results[0]["higher_parent_complete_replay_sha256"],
    "archive_sha256":results[0]["archive_sha256"],"workers":[{"parts":r["replayed_exception_root_parts"],
    "replayed":r["successful_original_trees_replayed"]}for r in results]},indent=2))
