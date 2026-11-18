"""Independently check the public receipt archive against all original v9 input trees."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile

parser=argparse.ArgumentParser(description=__doc__)
for key in ("checkout","certdir","table","archive","temp_parent"):
    parser.add_argument(key,type=Path)
parser.add_argument("--partial",action="store_true",help="allow partially completed demonstration archives")
args=parser.parse_args()
source=Path(__file__).resolve().parent
sha=lambda data:hashlib.sha256(data).hexdigest()
def digest_format(x):
    return isinstance(x,str) and len(x)==64 and all(c in "0123456789abcdef" for c in x)
with tempfile.TemporaryDirectory(dir=args.temp_parent) as scratch, ZipFile(args.archive) as z:
    out=Path(scratch)
    ledger=json.loads(z.read("coverage.json"))
    provenance=z.read("provenance.json")
    assert sha(provenance)==ledger["public_provenance_sha256"]
    origin=json.loads(provenance)
    assert origin["publication_script_sha256"]==sha((source/"publish_tomar_v3_archive.py").read_bytes())
    assert all(digest_format(origin[k]) for k in ("source_archive_sha256","source_coverage_sha256"))
    assert len(origin["rows"])==len(ledger["rows"])
    expected={"coverage.json","provenance.json"}
    hosts=set()
    for row,old in zip(ledger["rows"],origin["rows"]):
        p,host=row[:2]; assert host in ("local","macbook") and type(p) is int and 0<=p<2048
        assert old[0]==p and old[2]==row[-2] and all(digest_format(s) for s in old[1:])
        hosts.add(host)
        stem=f"part{p:04}"
        names=(f"{host}/integer/{stem}.json",f"{host}/integer/{stem}.npz",f"{host}/arb/{stem}.json")
        for name,d in zip(names,row[-3:]):
            assert name not in expected
            expected.add(name)
            b=z.read(name);assert sha(b)==d
            path=out/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b)
        native=json.loads(z.read(names[0]))
        assert [v for v in native["compiler_flags"] if v.startswith("-fmodules-cache-path=")]==["-fmodules-cache-path=<CLANG_MODULE_CACHE>"]
    assert set(z.namelist())==expected and len(z.namelist())==len(expected)
    for name in expected:
        if name.endswith(".json"):
            text=z.read(name)
            assert all(term not in text for term in (b"/Users/",b"/home/",b"/Volumes/"))
    cmd=[sys.executable,"-I","-B",str(source/"audit_tomar_v3_coverage.py"),args.checkout,args.certdir,args.table]
    for host in sorted(hosts):
        cmd.extend(["--pair",host,out/host/"integer",out/host/"arb"])
    if not args.partial:cmd.append("--require-complete")
    actual=json.loads(subprocess.check_output(cmd))
    # Complete published successes are self-contained. Transient native-only
    # parts may have been visible when a partial development archive was made.
    compare=("target","certified_global_11d_replay","certified_root_count","original_leaves","integer_closed",
             "arb_closed","arb_failed","missing_root_count","exterior_area_lower_bounds","native_identity","arb_identity","columns","rows")
    assert all(actual[k]==ledger[k] for k in compare)
    assert args.partial or (actual["certified_global_11d_replay"] and ledger["native_waiting_for_arb"]==[])
print(json.dumps({"all_roots":actual["certified_global_11d_replay"],"roots":actual["certified_root_count"],
    "original_leaves":actual["original_leaves"],"integer_closed":actual["integer_closed"],
    "arb_closed":actual["arb_closed"],"failed":actual["arb_failed"],"public_archive_sha256":sha(args.archive.read_bytes()),
    "private_archive_commitment":origin["source_archive_sha256"],"identity_path_scan":"passed"},indent=2))
