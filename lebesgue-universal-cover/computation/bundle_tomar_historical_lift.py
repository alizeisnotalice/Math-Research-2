"""Bundle a split statement with historical strong leaves and freshly checked weak leaves."""
import argparse,json,hashlib
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED
p=argparse.ArgumentParser(description=__doc__)
for name in ("selection","first","extra_selection","extra","output"):p.add_argument(name,type=Path)
a=p.parse_args();tmp=a.output.with_suffix(".partial")
files=[("selection.json",a.selection),("deeper-selection.json",a.extra_selection)]
for label,folder in (("first",a.first),("deeper",a.extra)):
 files += [(label+"/"+f.name,f) for f in sorted(folder.glob("part*.json"))]
with ZipFile(tmp,"w") as archive:
 for label,f in files:
  data=f.read_bytes();assert all(w not in data for w in (b"/Users/",b"/home/",b"/Volumes/"))
  ent=ZipInfo(label);ent.external_attr=0o100644<<16
  archive.writestr(ent,data,compress_type=ZIP_DEFLATED,compresslevel=9)
with ZipFile(tmp) as archive:assert archive.testzip() is None and len(archive.namelist())==len(files)
tmp.replace(a.output)
print(json.dumps({"bytes":a.output.stat().st_size,"members":len(files),"sha256":hashlib.sha256(a.output.read_bytes()).hexdigest()}))
