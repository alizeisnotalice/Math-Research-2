"""Seal the primary stronger-score partition, exact residual selection and all final checks."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

p=argparse.ArgumentParser(description=__doc__)
for k in ("primary","native","selection","extra","output"):p.add_argument(k,type=Path)
a=p.parse_args()
def sha(x):return hashlib.sha256(x).hexdigest()
report=json.loads(a.primary.read_bytes());sel=json.loads(a.selection.read_bytes())
assert sel["primary_sha256"]==sha(a.primary.read_bytes())
files={"primary.json":a.primary,"selection.json":a.selection}
for part,h1,h2 in report["higher_receipt_hashes"]:
    stem=f"part{part:04}"
    for suffix,h in ((".json",h1),(".npz",h2)):
        f=a.native/(stem+suffix);assert sha(f.read_bytes())==h;files["first/"+f.name]=f
for k,pos in sel["retained"].items():
    f=a.extra/f"part{int(k):04}.json";doc=json.loads(f.read_bytes())
    assert doc["selection_sha256"]==sha(a.selection.read_bytes()) and doc["tested"]==doc["closed"]==len(pos) and doc["remaining"]==[]
    files["extra/"+f.name]=f
tmp=a.output.with_suffix(".partial")
with ZipFile(tmp,"w") as output:
    for name,f in sorted(files.items()):
        data=f.read_bytes()
        if f.suffix==".json":assert all(s not in data for s in (b"/Users/",b"/home/",b"/Volumes/"))
        entry=ZipInfo(name);entry.external_attr=0o100644<<16
        output.writestr(entry,data,compress_type=ZIP_DEFLATED,compresslevel=9)
with ZipFile(tmp) as z:assert set(z.namelist())==set(files) and len(z.namelist())==len(files) and z.testzip() is None
assert tmp.stat().st_size<10*1024**2
tmp.replace(a.output)
print(json.dumps({"sha256":sha(a.output.read_bytes()),"members":len(files),"bytes":a.output.stat().st_size}))
