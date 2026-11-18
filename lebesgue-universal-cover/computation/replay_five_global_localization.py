"""Replay five-body high-threshold localization on original preorder partitions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import ctypes
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tomar_five_orbit_slice import CUT, leaf_mask
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table

THRESHOLD = Fraction("0.837295006222")
Q = 2**28
FINE = 2**44
STAGES = ((2, 127), (1, 1023), (1, 8191), (1, 65535))
TABLE_SHA = "fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
SOURCE = HERE / "five_grouped_localization.cpp"
INDEX = HERE / "five-orbit-slice.json"
FROZEN = HERE / "five-global-localization.json"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def make_table(path, cert):
    subprocess.run([sys.executable, "-I", "-B", HERE / "tomar_fan_table.py", cert, path,
                    "--subdiv", "32"], check=True)
    return path


def retained_envelope(codes, positions, lower, upper):
    """Obtain the outer coordinate box of the retained original closed leaves."""
    if not positions:
        return None
    end = [0] * len(codes)
    for i in range(len(codes) - 1, -1, -1):
        end[i] = i + 1 if codes[i] >= 253 else end[end[i + 1]]
    assert end[0] == len(codes)
    glo, ghi = upper.copy(), lower.copy()
    for pos in positions:
        lo, hi = lower.copy(), upper.copy()
        i = 0
        while i != pos:
            dim = int(codes[i])
            assert 0 <= dim < 11 and i < pos < end[i] and (lo[dim] + hi[dim]) % 2 == 0
            mid = (lo[dim] + hi[dim]) // 2
            if pos < end[i + 1]:
                hi[dim] = mid
                i += 1
            else:
                lo[dim] = mid
                i = end[i + 1]
        assert codes[i] >= 253
        glo, ghi = np.minimum(glo, lo), np.maximum(ghi, hi)
    return [glo.tolist(), ghi.tolist()]


def run():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cert", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("--table", type=Path, help="reuse an existing table and verify every point")
    p.add_argument("--parts", help="comma-separated root numbers for a partial replay")
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    assert __debug__ and 1 <= args.workers <= 32
    work = args.work.resolve()
    work.mkdir(exist_ok=False)
    cache = args.cache.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    cert = args.cert.resolve()
    assert digest(cert / "SHA256SUMS") == PINNED_MANIFEST
    manifest = {name.removeprefix("./"): val for val, name in
                (line.split() for line in (cert / "SHA256SUMS").read_text().splitlines())}
    assert digest(INDEX) == "ec483a5d64c86fc1d305c63081e08d66fd128bb2033622764ee26b2d45377e02"
    slice_index = json.loads(INDEX.read_bytes())
    assert Fraction(slice_index["cut"]) == CUT and slice_index["manifest"] == PINNED_MANIFEST
    assert digest(cert / "part0000.json") == manifest["part0000.json"]
    header = json.loads((cert / "part0000.json").read_bytes())
    table_path = args.table or make_table(work / "table.npz", cert)
    assert digest(table_path) == TABLE_SHA
    table, counts = verified_table(table_path, cert, header)
    assert counts.tolist() == [96, 96, 160, 224, 288]
    tables = {}
    for div, _ in STAGES:
        offsets = np.r_[0, np.cumsum(counts)]
        ids = [j for l, h in zip(offsets, offsets[1:]) for j in range(l, h, div)]
        sizes = np.asarray([(int(k)+div-1)//div for k in counts], np.int32)
        tables[div] = np.ascontiguousarray(table[:, ids, :]).ravel(), sizes
    radial_script = HERE / "certify_radial_affine_fans.py"
    radial = json.loads(subprocess.check_output([sys.executable, "-I", "-B", radial_script]))
    radii = [Fraction(r["radius"]) for r in radial["rows"]]
    assert Fraction(radial["target"]) == THRESHOLD
    assert radii == list(map(Fraction, [".06835", ".06650", ".06004", ".05699"]))
    with (work / "radial.json").open("w") as f:
        json.dump(radial, f)
    cpp = SOURCE.read_text()
    ceil = lambda v: -(-v.numerator // v.denominator)
    assert f"cut={ceil(CUT*FINE)};" in cpp
    assert ("r[4]={" + ",".join(str(ceil(r*FINE)) for r in radii) + "};") in cpp
    compiler = Path(os.environ.get("CXX", "/usr/bin/clang++"))
    flags = ["-O3", "-std=c++17", "-fsanitize=signed-integer-overflow",
             "-fsanitize-undefined-trap-on-error", "-fmodules-cache-path="+str(cache),
             "-shared", "-fPIC"]
    lib_path = work / "five-localization.dylib"
    subprocess.run([compiler, *flags, SOURCE, "-o", lib_path], check=True)
    fn = ctypes.CDLL(str(lib_path)).adaptive_trim
    v = lambda typ: np.ctypeslib.ndpointer(dtype=typ, ndim=1, flags="C_CONTIGUOUS")
    fn.argtypes = [v(np.uint8), ctypes.c_int, v(np.int64), v(np.int64), v(np.int64),
                   v(np.uint8), v(np.uint8), v(np.int64), v(np.int32), v(np.int32),
                   ctypes.c_int64, ctypes.c_int, ctypes.c_int]
    fn.restype = ctypes.c_int
    reference = json.loads(FROZEN.read_bytes())
    assert Fraction(reference["target"]) == THRESHOLD and reference["denominator"] == FINE
    assert len(reference["rows"]) == 2048
    exact_keys = {str(i) for i, row in enumerate(reference["rows"]) if row[-2]}
    assert set(reference["retained"]) == set(reference["envelopes"]) == exact_keys
    extremes = np.asarray(list(reference["envelopes"].values()), np.int64)
    assert reference["envelope"] == [extremes[:, 0].min(axis=0).tolist(), extremes[:, 1].max(axis=0).tolist()]
    from flint import arb
    assert arb.pi() / 15 < arb(CUT.numerator) / CUT.denominator
    for row, (k, radius) in enumerate(zip((3, 5, 7, 9), radii)):
        for d in (3 * row, 3 * row + 1):
            assert Fraction(header["lo"][d]) < -radius < radius < Fraction(header["hi"][d])
        if k in (7, 9):
            hi_angle = Fraction(header["hi"][3*row - 1])
            assert arb(hi_angle.numerator) / hi_angle.denominator > 2*arb.pi()/k
    assert Fraction(header["hi"][2]) > CUT and all(header["lo"][i] == 0 for i in (2, 5, 8))
    stages = reference["stages"]
    assert stages == [list(s) for s in STAGES]
    previous = reference["rows"]
    root_lo = np.array([int(Fraction(x) * FINE) for x in header["lo"]], np.int64)
    root_hi = np.array([int(Fraction(x) * FINE) for x in header["hi"]], np.int64)
    assert (root_lo+root_hi).tolist() == [2*int(v) for v in (root_lo+root_hi)//2]
    middle = (root_lo + root_hi)//2
    goal = int(2*Q*Q*THRESHOLD) + 1
    requested = list(range(2048)) if args.parts is None else sorted(set(map(int, args.parts.split(","))))
    assert requested and all(0 <= k < 2048 for k in requested)

    def process(part):
        stem = f"part{part:04}"
        src = cert / (stem + ".npy")
        meta_path = cert / (stem + ".json")
        assert digest(src) == manifest[src.name] and digest(meta_path) == manifest[meta_path.name]
        meta = json.loads(meta_path.read_bytes())
        assert meta["part"] == part and meta["split0"] == [2]*11 and meta["nparts"] == 2048
        assert meta["ns"] == [3, 5, 7, 9] and all(meta[k] == header[k] for k in ("lo","hi"))
        codes = np.load(src, allow_pickle=False)
        labels = leaf_mask(codes, part, Fraction(header["hi"][2]))
        srow = slice_index["rows"][part]
        assert srow[0] == part and sha256(labels.tobytes()).hexdigest() == srow[-1]
        assert np.bincount(labels, minlength=4)[1:].tolist() == srow[1:4]
        active = ((labels == 1) | (labels == 2)).astype(np.uint8)
        assert int(active.sum()) == srow[1] + srow[2]
        bits = np.array(np.unravel_index(part, (2,)*11), bool)
        lo, hi = root_lo.copy(), root_hi.copy()
        lo[bits], hi[~bits] = middle[bits], middle[~bits]
        row = [part, int(active.sum())]
        for div, budget in STAGES:
            new = np.zeros_like(active)
            stats = np.zeros(6, np.int64)
            if active.any():
                values, sizes = tables[div]
                assert fn(codes, len(codes), lo, hi, root_hi, active, new, stats, values,
                          sizes, goal, 48, budget) == 0
            assert np.all(new <= active) and new.sum()+stats[2] == active.sum()
            active = new
            row.extend([int(active.sum()), int(stats[0])])
        assert previous[part] == row
        final = np.flatnonzero(active).tolist()
        assert final == reference["retained"].get(str(part), [])
        assert retained_envelope(codes, final, lo, hi) == reference["envelopes"].get(str(part))
        return row

    ordered = sorted(requested, key=lambda p: -sum(slice_index["rows"][p][1:3]))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = sorted(pool.map(process, ordered))
    total = np.sum(np.array(rows, dtype=np.int64)[:, 1:], axis=0).tolist()
    if args.parts is None:
        assert total == reference["totals"]
    frozen = dict(complete=args.parts is None, parts=len(rows), totals=total, goal=goal,
                  threshold=str(THRESHOLD), python=platform.python_version(),
                  compiler_version=subprocess.check_output([compiler,"--version"],text=True).splitlines()[0],
                  compiler_flags=flags, library_sha256=digest(lib_path), source_sha256=digest(SOURCE),
                  driver_sha256=digest(Path(__file__)), expected_sha256=digest(FROZEN),
                  certified_table_sha256=digest(table_path), radial_verification_sha256=digest(radial_script),
                  source_manifest=PINNED_MANIFEST, slice_index_sha256=digest(INDEX))
    (work / "replay.json").write_text(json.dumps(frozen, indent=2)+"\n")
    print(json.dumps(frozen))


if __name__ == "__main__":
    run()
