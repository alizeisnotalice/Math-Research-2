"""Audit complete original-part receipts and exact residual membership for a V3 replay."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("checkout", type=Path)
parser.add_argument("certdir", type=Path)
parser.add_argument("table", type=Path)
parser.add_argument("--pair", nargs=3, action="append", required=True, metavar=("HOST", "NATIVE", "ARB"))
parser.add_argument("--require-complete", action="store_true")
args = parser.parse_args()
assert __debug__, "assertions must be enabled"
ctx.prec = 128

def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def read(path):
    return json.loads(path.read_bytes())

target = "0.83661"
table_digest = "0c5bcf5dc69d8e8398dfbdb209f1563d5623523ffd5bca21ae6da58b6d92d88c"
library_digest = "ae3814e8832165e2ead77af3def4f03948dddf78603f06f51ac8b23cbe7be525"
manifest = args.certdir / "SHA256SUMS"
assert sha(manifest) == PINNED_MANIFEST
pinned = {name.removeprefix("./"): value for value, name in
          (line.split() for line in manifest.read_text().splitlines())}
first_path = args.certdir / "part0000.json"
assert sha(first_path) == pinned[first_path.name]
first = read(first_path)
assert first["ns"] == [3, 5, 7, 9] and first["split0"] == [2]*11
assert len(first["lo"]) == len(first["hi"]) == 11 and first["nparts"] == 2048
assert sha(args.table) == table_digest
verified_table(args.table, args.certdir, first)
upstream = args.checkout / "checker" / "lbcheck.py"
revision = "10af98606c98ca9ffd0d9f402cae1e2efe74b792"
assert subprocess.check_output(["git", "-C", str(args.checkout), "show", revision+":checker/lbcheck.py"]) == upstream.read_bytes()
native_binding = {
    "cpp_source_sha256": sha(HERE / "tomar_integer_fan.cpp"),
    "generator_sha256": sha(HERE / "tomar_fan_table.py"),
    "verifier_sha256": sha(HERE / "verify_tomar_fan_table.py"),
    "table_sha256": table_digest, "library_sha256": library_digest,
    "input_manifest_sha256": PINNED_MANIFEST, "prior_unresolved_sha256": None,
    "initial_mask": "all original nodes with code 253, 254 or 255",
    "target": target, "scale": 2**28, "table_grid": 1024
}
arb_binding = {
    "target": target, "source_sha256": sha(HERE / "finish_tomar_v3_parts.py"),
    "decoder_sha256": sha(HERE / "tomar_leaf_box.py"), "upstream_checker_sha256": sha(upstream)
}
outside = []
for i, n in enumerate(first["ns"]):
    x, y = 3*i, 3*i+1
    d = min(-first["lo"][x], first["hi"][x], -first["lo"][y], first["hi"][y])
    assert d > 0
    radius = 1 - 1/(2*(arb.pi()/(2*n)).cos())
    exterior = arb.pi()*(arb(1)/4 + radius**2)/2 + arb(d)*(arb(1)/2+radius)
    assert exterior > arb(target)
    outside.append(str(exterior))
# Fix the B3 axis by global rotation, then reflect the pentagon if necessary.
# The remaining B7 and B9 orientations use the full period.
for coordinate, symmetry, numerator in [(2, 5, 1), (5, 7, 2), (8, 9, 2)]:
    assert first["lo"][coordinate] == 0 and arb(first["hi"][coordinate]) > numerator*arb.pi()/symmetry

rows = {}
pending = []
for host, path_native, path_arb in args.pair:
    assert host in ("local", "macbook")
    driver = "run_tomar_integer_parts.py" if host == "local" else "run_tomar_portable_parts.py"
    driver_digest = sha(HERE / driver)
    native_dir, arb_dir = Path(path_native), Path(path_arb)
    assert native_dir.is_dir() and arb_dir.is_dir()
    for stage in sorted(native_dir.glob("part*.json")):
        stem = stage.stem
        part = int(stem.removeprefix("part"))
        assert 0 <= part < 2048 and stem == f"part{part:04}" and part not in rows and part not in pending
        end = arb_dir / (stem+".json")
        if not end.is_file():
            pending.append(part)
            continue
        meta_path, code_path = args.certdir/(stem+".json"), args.certdir/(stem+".npy")
        native_bytes, arb_bytes = stage.read_bytes(), end.read_bytes()
        native_sha = hashlib.sha256(native_bytes).hexdigest()
        arb_sha = hashlib.sha256(arb_bytes).hexdigest()
        meta, record, finish = read(meta_path), json.loads(native_bytes), json.loads(arb_bytes)
        assert meta["part"] == record["part"] == finish["part"] == part
        assert all(meta[k] == first[k] for k in ("lo","hi","ns","split0","nparts","T","eps","N"))
        assert sha(meta_path) == pinned[meta_path.name] == record["header_sha256"]
        assert sha(code_path) == pinned[code_path.name] == record["codes_sha256"]
        assert all(record[k] == v for k, v in native_binding.items()) and record["driver_sha256"] == driver_digest
        assert all(finish[k] == v for k, v in arb_binding.items())
        assert native_sha == finish["input_receipt_sha256"]
        masks = native_dir / (stem+".npz")
        assert sha(masks) == record["remaining_sha256"] == finish["prior_unresolved_sha256"]
        codes = np.load(code_path, allow_pickle=False)
        assert codes.dtype == np.uint8 and codes.ndim == 1
        assert np.all((codes < 11) | (codes >= 253))
        balance = 1 + np.cumsum(np.where(codes < 11, 1, -1))
        assert np.all(balance[:-1] > 0) and balance[-1] == 0
        original = int(np.count_nonzero(codes >= 253))
        assert original == record["all_original"] == record["tested"] == finish["all_original"]
        assert original == record["closed"]+record["remaining"] and record["closed"] == finish["closed_by_stage"]
        with np.load(masks, allow_pickle=False) as bundle:
            assert bundle.files == ["pos"]
            pos = bundle["pos"]
        assert pos.dtype == np.int32 and pos.ndim == 1 and np.array_equal(pos, np.unique(pos))
        assert np.all((pos >= 0) & (pos < len(codes))) and np.all(codes[pos] >= 253)
        assert len(pos) == record["remaining"] == finish["arb_attempts"] == len(finish["leaves"])
        assert record["remaining_codes"] == {str(k):int(np.count_nonzero(codes[pos] == k)) for k in (253,254,255)}
        failures = []
        for index, node in zip(pos.tolist(), finish["leaves"]):
            index2, tag, calls, depth, passed = node
            assert index2 == index and tag == int(codes[index]) and type(passed) is bool
            assert 0 <= depth <= finish["depth"] and calls >= 1+depth
            if not passed:
                failures.append(index)
        assert failures == finish["unresolved"] and finish["arb_closed"] == len(pos)-len(failures)
        rows[part] = [host, original, record["closed"], finish["arb_closed"], len(failures),
                      native_sha, sha(masks), arb_sha]
parts = sorted(rows)
original, native, checked, failed = (sum(rows[i][j] for i in parts) for j in (1,2,3,4))
full = parts == list(range(2048)) and original == 16193402 and failed == 0 and not pending
assert not args.require_complete or full, "original roots or positive area checks are missing"
print(json.dumps({"target":target, "certified_global_11d_replay":full, "certified_root_count":sum(rows[i][4]==0 for i in parts),
    "original_leaves":original, "integer_closed":native, "arb_closed":checked, "arb_failed":failed,
    "native_waiting_for_arb":sorted(pending), "missing_root_count":2048-len(parts), "exterior_area_lower_bounds":outside,
    "source_sha256":sha(Path(__file__)), "native_identity":native_binding, "arb_identity":arb_binding,
    "columns":["part","host","original","integer_closed","arb_closed","failed","native_sha","mask_sha","arb_sha"],
    "rows":[[i,*rows[i]] for i in parts]},indent=2))
