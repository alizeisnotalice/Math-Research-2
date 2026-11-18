"""Seal every fixed-fan and dynamic-affine receipt of the 0.83725 localization lift."""
import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("archive",type=Path)
p.add_argument("--level",action="append",nargs=2,metavar=("SELECTION","RECEIPTS"),required=True)
a=p.parse_args();assert len(a.level)==5
tmp=a.archive.with_suffix(".partial");members=0
with ZipFile(tmp,"w")as z:
    for stage,(selection,receipts)in enumerate(a.level,1):
        selection=Path(selection);receipts=Path(receipts);s=json.loads(selection.read_bytes())
        files={"selection.json":selection}
        for k in sorted(map(int,s["retained"])):
            name=f"part{k:04}.json";f=receipts/name;d=json.loads(f.read_bytes())
            assert d["part"]==k and d["selection_sha256"]==hashlib.sha256(selection.read_bytes()).hexdigest()
            if stage==5:assert all(not x["remaining"]for x in d["rows"])
            files[name]=f
        for name,f in files.items():
            value=f.read_bytes();assert all(x not in value for x in (b"/Users/",b"/Volumes/",b"/home/"))
            ent=ZipInfo(f"level{stage}/{name}");ent.external_attr=0o100644<<16
            z.writestr(ent,value,compress_type=ZIP_DEFLATED,compresslevel=9);members+=1
with ZipFile(tmp)as z:assert z.testzip() is None and len(z.namelist())==members
tmp.replace(a.archive);print(json.dumps({"members":members,"bytes":a.archive.stat().st_size,
    "sha256":hashlib.sha256(a.archive.read_bytes()).hexdigest()}))
