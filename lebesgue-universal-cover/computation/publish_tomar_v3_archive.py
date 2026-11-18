"""Produce public receipts while retaining commitments to every private source receipt."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("private", type=Path)
parser.add_argument("public", type=Path)
args = parser.parse_args()
assert args.private.resolve() != args.public.resolve()

def digest(data):
    return hashlib.sha256(data).hexdigest()

def encode(obj):
    return (json.dumps(obj,indent=2)+"\n").encode()

def add(z, name, data):
    if name.endswith(".json"):
        # Absolute OS paths would reproduce host account or storage locations.
        assert b"/Users/" not in data and b"/home/" not in data and b"/Volumes/" not in data
    entry=ZipInfo(name)
    entry.external_attr=0o100644 << 16
    z.writestr(entry,data,compress_type=ZIP_DEFLATED,compresslevel=9)

tmp = args.public.with_suffix(".partial")
with ZipFile(args.private) as src, ZipFile(tmp,"w") as dst:
    private_coverage=src.read("coverage.json")
    private=json.loads(private_coverage)
    public=copy.deepcopy(private)
    original_rows=[]
    expected={"coverage.json"}
    for row in public["rows"]:
        part,host = row[:2]
        stem=f"part{part:04}"
        n,m,a=(f"{host}/integer/{stem}.json",f"{host}/integer/{stem}.npz",f"{host}/arb/{stem}.json")
        expected.update([n,m,a])
        raw_native,mask,raw_arb=map(src.read,(n,m,a))
        assert [digest(x) for x in (raw_native,mask,raw_arb)]==row[-3:]
        native,arb=json.loads(raw_native),json.loads(raw_arb)
        assert arb["input_receipt_sha256"]==row[-3]
        assert native["remaining_sha256"]==arb["prior_unresolved_sha256"]==row[-2]
        flags=native["compiler_flags"]
        prefix="-fmodules-cache-path="
        match=[i for i,s in enumerate(flags) if s.startswith(prefix)]
        assert len(match)==1 and flags[match[0]][len(prefix):].startswith("/")
        flags[match[0]]=prefix+"<CLANG_MODULE_CACHE>"
        new_native=encode(native)
        arb["input_receipt_sha256"]=digest(new_native)
        new_arb=encode(arb)
        original_rows.append([part,*row[-3:]])
        row[-3:]=[digest(new_native),digest(mask),digest(new_arb)]
        for name,data in ((n,new_native),(m,mask),(a,new_arb)):
            add(dst,name,data)
    assert len(src.namelist())==len(expected) and set(src.namelist())==expected
    provenance={"transformation":"Only the module cache directory in compiler_flags is replaced; the Arb receipt and coverage ledger bind the resulting bytes.",
        "source_archive_sha256":digest(args.private.read_bytes()),"source_coverage_sha256":digest(private_coverage),
        "publication_script_sha256":digest(Path(__file__).read_bytes()),
        "columns":["part","source_native_sha256","source_mask_sha256","source_arb_sha256"],"rows":original_rows}
    source_bytes=encode(provenance)
    public["public_provenance_sha256"]=digest(source_bytes)
    add(dst,"provenance.json",source_bytes)
    add(dst,"coverage.json",encode(public))
with ZipFile(tmp) as z:
    assert z.testzip() is None
    r=json.loads(z.read("coverage.json"))
    assert digest(z.read("provenance.json"))==r["public_provenance_sha256"]
    assert len(z.namelist())==3*len(r["rows"])+2
tmp.replace(args.public)
print(json.dumps({"roots":len(r["rows"]),"public_bytes":args.public.stat().st_size,
                  "private_sha256":provenance["source_archive_sha256"],
                  "public_sha256":digest(args.public.read_bytes())},indent=2))
