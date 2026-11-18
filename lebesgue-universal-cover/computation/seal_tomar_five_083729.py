"""Package every candidate-placement receipt for the 501 high-localization exceptions at 0.83729."""
import argparse,hashlib,json
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED
p=argparse.ArgumentParser(description=__doc__)
for k in ("high","active","first","archive"):p.add_argument(k,type=Path)
p.add_argument("--five",type=Path,action="append",required=True)
a=p.parse_args()
def sha(b):return hashlib.sha256(b).hexdigest()
active=json.loads(a.active.read_bytes());high=json.loads(a.high.read_bytes())
assert sha(a.high.read_bytes())==active["high_selection_sha256"]=="d85941b19b47cda83fc86268daa50569db2d80118b7fc722d4bf90adfa457504"
files={"active.json":a.active};done=set()
for part,pos,h in active["first_receipts"]:
    name=f"part{part:04}-pos{pos}.json";f=a.first/name;r=json.loads(f.read_bytes())
    assert sha(f.read_bytes())==h and r["part"]==part and r["row"]["label"]==pos
    files["first/"+name]=f
    if not r["row"]["remaining"]:done.add((part,pos))
pool={f.name:f for folder in a.five for f in folder.glob("part*-pos*.json")}
assert len(pool)==sum(len(list(folder.glob("part*-pos*.json")))for folder in a.five)
for k,ids in active["retained"].items():
    for pos in ids:
        n=f"part{int(k):04}-pos{pos}.json";f=pool.pop(n);r=json.loads(f.read_bytes())
        assert r["part"]==int(k) and r["row"]["label"]==pos and not r["row"]["remaining"]
        assert r["active_sha256"]==sha(a.active.read_bytes())
        files["five/"+n]=f;done.add((int(k),pos))
assert not pool and len(done)==501 and done=={(int(k),i)for k,v in high["retained"].items()for i in v}
tmp=a.archive.with_suffix(".partial")
with ZipFile(tmp,"w")as z:
    for name,f in sorted(files.items()):
        raw=f.read_bytes();assert all(b not in raw for b in (b"/Users/",b"/home/",b"/Volumes/"))
        e=ZipInfo(name);e.external_attr=0o100644<<16;z.writestr(e,raw,compress_type=ZIP_DEFLATED,compresslevel=9)
with ZipFile(tmp)as z:assert z.testzip() is None and set(z.namelist())==set(files)
assert tmp.stat().st_size<10*1024**2
tmp.replace(a.archive);print(json.dumps({"sha256":sha(a.archive.read_bytes()),"members":len(files),"bytes":a.archive.stat().st_size}))
