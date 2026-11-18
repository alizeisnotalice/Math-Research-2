"""Reconstruct the D3 slice from the original trees and audit both grouped-fan receipt stages."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_tomar_grouped_residues import TARGET, TABLE, DEPTH, CALLS, sha
from tomar_d3_mask import leaf_mask
from verify_tomar_fan_table import verified_table, PINNED_MANIFEST

INDEX = "ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
PRESET = ["-O3", "-ftrapv", "-std=c++17"]
LIBS = [("tomar_grouped_scored_fan.cpp", "libtomar_grouped_scored.dylib"),
        ("tomar_grouped_fan.cpp", "libtomar_grouped_residue.dylib")]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("certdir", "table", "slice_index", "native", "supplement", "scratch", "cache"):
        p.add_argument(key, type=Path)
    p.add_argument("--partial", action="store_true")
    args = p.parse_args()
    assert sha(args.table) == TABLE and sha(args.slice_index) == INDEX
    manifest = args.certdir / "SHA256SUMS"
    assert sha(manifest) == PINNED_MANIFEST
    pinned = {n.removeprefix("./"): h for h, n in (r.split() for r in manifest.read_text().splitlines())}
    meta = json.loads((args.certdir / "part0000.json").read_bytes())
    verified_table(args.table, args.certdir, meta)
    rows = json.loads(args.slice_index.read_bytes())["rows"]
    assert len(rows) == 2048 and [r[0] for r in rows] == list(range(2048))
    required = {r[0] for r in rows if r[1] + r[2] > 0}
    assert len(required) == 1024 and sum(r[1] + r[2] for r in rows) == 6710314
    assert sum(r[3] for r in rows) == 9483088
    parts = {int(f.stem[4:]) for f in args.supplement.glob("part*.json")}
    native_parts = {int(f.stem[4:]) for f in args.native.glob("part*.json")}
    assert parts and parts <= native_parts <= required
    if not args.partial:
        assert parts == native_parts == required
    with tempfile.TemporaryDirectory(dir=args.scratch) as tmp:
        libs = []
        for source, name in LIBS:
            output = Path(tmp) / name
            subprocess.run(["/usr/bin/clang++", *PRESET, f"-fmodules-cache-path={args.cache}", "-shared",
                            "-fPIC", f"-Wl,-install_name,@rpath/{name}", HERE / source, "-o", output], check=True)
            libs.append(sha(output))
    compiler = Path("/usr/bin/clang++")
    flags = [*PRESET, "-fmodules-cache-path=<CLANG_MODULE_CACHE>", "-shared", "-fPIC",
             "-Wl,-install_name,@rpath/libtomar_grouped_scored.dylib"]
    expected = {"target": TARGET, "scale": 2**28, "table_grid": 1024, "table_sha256": TABLE,
                "slice_index_sha256": INDEX, "slice_selector_sha256": sha(HERE / "tomar_d3_mask.py"),
                "max_depth": 30, "max_calls": 4095, "library_sha256": libs[0], "compiler_sha256": sha(compiler),
                "compiler_version": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
                "compiler_flags": flags, "input_manifest_sha256": PINNED_MANIFEST,
                "generator_sha256": sha(HERE / "tomar_fan_table.py"), "verifier_sha256": sha(HERE / "verify_tomar_fan_table.py"),
                "driver_sha256": sha(HERE / "run_tomar_grouped_slice.py"), "cpp_source_sha256": sha(HERE / LIBS[0][0]),
                "prior_unresolved_sha256": None, "initial_mask": "preorder leaf_mask(codes,part,original_angle_high) in {1,2}"}
    expected_extra = {"target": TARGET, "table_sha256": TABLE, "driver_sha256": sha(HERE / "check_tomar_grouped_residues.py"),
                      "source_sha256": sha(HERE / LIBS[1][0]), "decoder_sha256": sha(HERE / "tomar_leaf_box.py"),
                      "library_sha256": libs[1], "depth": DEPTH, "calls": CALLS}
    totals = Counter()
    levels = Counter()
    committed = []
    scan = sorted(parts) if args.partial else range(2048)
    for part in scan:
        stem = f"part{part:04}"
        header = args.certdir / (stem + ".json")
        tree = args.certdir / (stem + ".npy")
        assert sha(header) == pinned[header.name] and sha(tree) == pinned[tree.name]
        h = json.loads(header.read_bytes())
        assert h["part"] == part and h["ns"] == [3, 5, 7, 9] and h["split0"] == [2]*11 and h["nparts"] == 2048
        assert h["lo"] == meta["lo"] and h["hi"] == meta["hi"]
        codes = np.load(tree, allow_pickle=False)
        mask = leaf_mask(codes, part, Fraction(h["hi"][2]))
        row = rows[part]
        assert [np.count_nonzero(mask == t) for t in (1, 2, 3)] == row[1:4]
        assert hashlib.sha256(mask.tobytes()).hexdigest() == row[-1]
        totals["raw_leaves_scanned"] += int(np.count_nonzero(codes >= 253))
        if part not in parts:
            assert row[1] + row[2] == 0
            continue
        primary = args.native / header.name
        npz = primary.with_suffix(".npz")
        deeper = args.supplement / header.name
        first, extra = [json.loads(f.read_bytes()) for f in (primary, deeper)]
        for key, val in expected.items():
            assert first[key] == val, (part, key)
        for key, val in expected_extra.items():
            assert extra[key] == val, (part, key)
        for receipt in [first, extra]:
            assert receipt["part"] == part and receipt["header_sha256"] == pinned[header.name] and receipt["codes_sha256"] == pinned[tree.name]
        assert extra["parent_sha256"] == sha(primary) and extra["positions_sha256"] == first["remaining_sha256"] == sha(npz)
        assert first["slice_mask_sha256"] == row[-1] and first["full_raw_leaf_count"] == sum(row[1:4])
        ids = np.flatnonzero((mask == 1) | (mask == 2))
        with np.load(npz, allow_pickle=False) as saved:
            assert set(saved.files) == {"pos", "scores"}
            pos, score = saved["pos"], saved["scores"]
        assert score.dtype == np.uint16 and score.shape == ids.shape and pos.dtype == np.int32
        assert np.array_equal(ids[score == 0], pos)
        assert first["tested"] == first["all_original"] == len(ids) and first["closed"] == np.count_nonzero(score)
        assert first["remaining"] == len(pos) == extra["checked"]
        assert first["remaining_codes"] == {str(k): int(np.count_nonzero(codes[pos] == k)) for k in (253, 254, 255)}
        counts = {str(k): int(np.count_nonzero(score >= k - 8371000 + 1)) for k in (8371000, 8371500, 8372000, 8372500, 8372700, 8372900)}
        assert counts == first["scores_at_1e7_targets"]
        levels.update(counts)
        entries = extra["leaves"]
        assert [r[:2] for r in entries] == [[int(i), int(codes[i])] for i in pos]
        assert extra["closed"] == sum(r[-1] for r in entries) == len(entries)
        for _, _, calls, depth, terminal, ok in entries:
            assert ok == 1 and 1 <= terminal <= calls <= CALLS and 0 <= depth <= DEPTH
        totals.update({"canonical_leaves": len(ids), "native_closed": first["closed"], "supplement_closed": len(entries),
                       "native_calls": first["checked_boxes"], "supplement_calls": sum(r[2] for r in entries)})
        committed.append([part, sha(primary), sha(npz), sha(deeper)])
    assert totals["canonical_leaves"] == totals["native_closed"] + totals["supplement_closed"]
    if not args.partial:
        assert totals["canonical_leaves"] == 6710314 and totals["raw_leaves_scanned"] == 16193402
    ctx.prec = 128
    outside = []
    for i, n in enumerate(meta["ns"]):
        j = 3*i
        r = 1 - 1/(2*(arb.pi()/(2*n)).cos())
        dist = min(-meta["lo"][j], meta["hi"][j], -meta["lo"][j+1], meta["hi"][j+1])
        area = arb.pi()*(arb(1)/4 + r*r)/2 + arb(dist)*(arb(1)/2 + r)
        assert area > arb(TARGET)
        outside.append(str(area))
    print(json.dumps({"target": TARGET, "complete": parts == required and not args.partial, "parts": len(parts),
                      "input_tree_sha256": PINNED_MANIFEST, "slice_index_sha256": INDEX, "table_sha256": TABLE,
                      "historical_area_inputs": [], "totals": dict(totals), "stronger_native_leaf_counts": dict(levels),
                      "root_exterior": outside, "compiled_libraries": libs, "receipt_hashes": committed},
                     separators=(",", ":")))

if __name__ == "__main__":
    main()
