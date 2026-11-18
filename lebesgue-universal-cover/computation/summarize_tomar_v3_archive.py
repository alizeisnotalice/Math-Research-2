"""Summarize a final frozen proof archive and check the hashes of every archived member."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys
from zipfile import ZipFile

path = Path(sys.argv[1])
by_host = defaultdict(Counter)
codes = Counter()
depths = Counter()
compiler = set()
fallback_max_depth = 0
with ZipFile(path) as z:
    ledger = json.loads(z.read("coverage.json"))
    assert ledger["certified_global_11d_replay"]
    assert [row[0] for row in ledger["rows"]] == list(range(2048))
    expected_names = {"coverage.json"}
    if "public_provenance_sha256" in ledger:
        expected_names.add("provenance.json")
        assert hashlib.sha256(z.read("provenance.json")).hexdigest()==ledger["public_provenance_sha256"]
    for part, host, original, integer, arb, failed, native_sha, mask_sha, arb_sha in ledger["rows"]:
        assert original == integer+arb and failed == 0
        prefix = f"part{part:04}"
        names=(f"{host}/integer/{prefix}.json",f"{host}/integer/{prefix}.npz",f"{host}/arb/{prefix}.json")
        expected_names.update(names)
        files=[z.read(name) for name in names]
        assert [hashlib.sha256(data).hexdigest() for data in files] == [native_sha,mask_sha,arb_sha]
        n,a = json.loads(files[0]),json.loads(files[2])
        assert n["tested"]==original and a["arb_closed"]==arb and a["unresolved"]==[]
        assert a["input_receipt_sha256"]==native_sha and n["remaining_sha256"]==a["prior_unresolved_sha256"]==mask_sha
        b=by_host[host]
        b.update(roots=1,original=original,integer=integer,arb=arb,integer_box_calls=n["checked_boxes"],
                 integer_part_wall_milliseconds=round(n["wall_seconds"]*1000))
        compiler.add((n["compiler_version"],n["library_sha256"]))
        for _,tag,calls,depth,passed in a["leaves"]:
            assert passed
            b["arb_box_calls"]+=calls
            codes[tag]+=1
            depths[depth]+=1
            fallback_max_depth=max(fallback_max_depth,depth)
    assert set(z.namelist())==expected_names and len(z.namelist())==len(expected_names)
assert sum(b["original"] for b in by_host.values())==16193402
assert sum(b["arb"] for b in by_host.values())==sum(codes.values())==ledger["arb_closed"]
print(json.dumps({"target":ledger["target"],"global":True,"parts":2048,"original":ledger["original_leaves"],
    "integer_closed":ledger["integer_closed"],"arb_closed":ledger["arb_closed"],"remaining":0,
    "by_host":by_host,"fallback_source_codes":codes,"fallback_leaf_maximum_depth_distribution":depths,
    "fallback_maximum_depth":fallback_max_depth,"compiler_and_library":sorted(compiler),
    "root_exterior_area_lower_bounds":ledger["exterior_area_lower_bounds"],
    "archived_files":len(expected_names),"archive_bytes":path.stat().st_size,
    "archive_sha256":hashlib.sha256(path.read_bytes()).hexdigest()},indent=2))
