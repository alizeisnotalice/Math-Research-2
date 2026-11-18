"""Public audit of the 0.8372 five-body bound from its 0.8371 score archive and the final supplement."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

p=argparse.ArgumentParser(description=__doc__)
for key in ("cert","table","slice_index","parent_ledger","parent_index","archive","scratch","cache"):
    p.add_argument(key,type=Path)
p.add_argument("--replay-exceptions",action="store_true")
a=p.parse_args()
here=Path(__file__).resolve().parent
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def run(script,*args):return subprocess.check_output([sys.executable,"-I","-B",here/script,*args])
base=json.loads(run("audit_tomar_grouped_public.py",a.cert,a.table,a.slice_index,
                    a.parent_ledger,a.parent_index,a.scratch,a.cache))
assert base["complete"] and base["historical_area_inputs"]==[] and base["parts"]==1024
parent=json.loads(a.parent_ledger.read_bytes());parent_index=json.loads(a.parent_index.read_bytes())
assert parent["totals"]["canonical_leaves"]==6710314
with tempfile.TemporaryDirectory(dir=a.scratch) as tmp:
    d=Path(tmp)
    for k in ("parents","first","extra"): (d/k).mkdir()
    for archive in parent_index["archives"]:
        with ZipFile(a.parent_index.parent/archive["name"]) as z:
            for n in z.namelist():
                if n.startswith("native/"):
                    (d/"parents"/Path(n).name).write_bytes(z.read(n))
    names={"primary.json","selection.json"}
    with ZipFile(a.archive) as z:
        assert z.testzip() is None
        for f in ("primary.json","selection.json"):(d/f).write_bytes(z.read(f))
        stored=json.loads((d/"primary.json").read_bytes());sel=json.loads((d/"selection.json").read_bytes())
        assert stored["parts"]==1024 and stored["authenticated_parent_ledger_sha256"]==sha(a.parent_ledger)
        assert sel["primary_sha256"]==sha(d/"primary.json") and sel["parent_sha256"]==sha(a.parent_ledger)
        for k,h1,h2 in stored["higher_receipt_hashes"]:
            stem=f"part{k:04}"
            for suffix,h in ((".json",h1),(".npz",h2)):
                n="first/"+stem+suffix;names.add(n);data=z.read(n)
                assert hashlib.sha256(data).hexdigest()==h;(d/n).write_bytes(data)
        for k in sel["retained"]:
            n=f"extra/part{int(k):04}.json";names.add(n);(d/n).write_bytes(z.read(n))
        assert set(z.namelist())==names and len(z.namelist())==len(names)
    rebuilt=json.loads(run("audit_tomar_score_lift.py",a.parent_ledger,d/"parents",a.cert,d/"first",d,a.cache,
                           "--target","0.8372","--partial"))
    assert rebuilt==stored
    run("select_tomar_score_exceptions.py",a.parent_ledger,d/"primary.json",d/"first",d/"regenerated.json")
    assert (d/"regenerated.json").read_bytes()==(d/"selection.json").read_bytes()
    source=here/"tomar_fine_angle_fan.cpp";driver=here/"check_tomar_score_exceptions.py"
    binary=d/"libtomar_score_exceptions.dylib"
    subprocess.run(["/usr/bin/clang++","-O3","-ftrapv","-std=c++17",f"-fmodules-cache-path={a.cache}","-shared",
                    "-fPIC","-Wl,-install_name,@rpath/libtomar_score_exceptions.dylib",source,"-o",binary],check=True)
    if a.replay_exceptions:
        run("check_tomar_score_exceptions.py",d/"selection.json",a.cert,a.table,a.slice_index,source,d/"repeat",a.cache)
    extra=visits=0;depth=0
    for k,pos in sel["retained"].items():
        file=d/"extra"/f"part{int(k):04}.json";r=json.loads(file.read_bytes())
        for key,value in {"part":int(k),"target":"0.8372","input_kind":sel["kind"],"selection_sha256":sha(d/"selection.json"),
            "slice_index_sha256":sha(a.slice_index),"header_sha256":sha(a.cert/file.name),
            "codes_sha256":sha(a.cert/file.with_suffix(".npy").name),"source_sha256":sha(source),
            "library_sha256":sha(binary),"driver_sha256":sha(driver),"table_sha256":sha(a.table),
            "depth":100,"calls":16777216,"tested":len(pos),"closed":len(pos),"remaining":[]}.items():
            assert r[key]==value,(k,key)
        if a.replay_exceptions:assert file.read_bytes()==(d/"repeat"/file.name).read_bytes()
        assert 0<=r["maximum_depth"]<=100
        extra+=len(pos);visits+=r["checked_boxes"];depth=max(depth,r["maximum_depth"])
    totals=stored["totals"]
    assert extra==totals["remaining"] and totals["reused"]+totals["fresh_closed"]+extra==totals["canonical"]==6710314
print(json.dumps({"target":"0.8372","complete":True,"parts":1024,"canonical_leaves":6710314,
    "full_original_leaves":parent["totals"]["raw_leaves_scanned"],"from_stored_area_scores":totals["reused"],
    "fresh_primary":totals["fresh_closed"],"fresh_exceptions":extra,"primary_box_calls":totals["checked_boxes"],
    "exception_box_calls":visits,"exception_maximum_depth":depth,"exceptions_replayed":a.replay_exceptions,
    "parent_ledger_sha256":sha(a.parent_ledger),"parent_bundle_index_sha256":sha(a.parent_index),
    "archive_sha256":sha(a.archive),"archive_members":len(names),"archive_bytes":a.archive.stat().st_size,
    "historical_area_inputs":[]},indent=2))
