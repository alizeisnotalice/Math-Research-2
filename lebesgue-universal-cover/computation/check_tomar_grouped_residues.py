"""Check every outstanding position of the native five-body slice by a deeper binary partition."""
import argparse
import ctypes
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import verified_table, PINNED_MANIFEST

TARGET = "0.8371"
TABLE = "fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c"
FINE = 2**44
SCALE = 2**28
DEPTH = 60
CALLS = 1048575

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def scaled(xs):
    values = [FINE * Fraction(x) for x in xs]
    assert all(x.denominator == 1 for x in values)
    return np.array([int(x) for x in values], dtype=np.int64)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("certdir", "table", "native", "out", "cache"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--parts")
    args = parser.parse_args()
    assert sha(args.table) == TABLE
    manifest = args.certdir / "SHA256SUMS"
    assert sha(manifest) == PINNED_MANIFEST
    pinned = {name.removeprefix("./"): h for h, name in (r.split() for r in manifest.read_text().splitlines())}
    header = json.loads((args.certdir / "part0000.json").read_bytes())
    points, counts = verified_table(args.table, args.certdir, header)
    base = scaled(header["hi"])
    source = HERE / "tomar_grouped_fan.cpp"
    args.out.mkdir(parents=True, exist_ok=True)
    keys = sorted(args.native.glob("part*.json"))
    if args.parts:
        subset = set(map(int, args.parts.split(",")))
        keys = [file for file in keys if int(file.stem[4:]) in subset]
        assert len(keys) == len(subset)
    with tempfile.TemporaryDirectory(dir=args.out) as tmp:
        binary = Path(tmp) / "libtomar_grouped_residue.dylib"
        subprocess.run(["/usr/bin/clang++", "-O3", "-ftrapv", "-std=c++17",
                        f"-fmodules-cache-path={args.cache}", "-shared", "-fPIC",
                        "-Wl,-install_name,@rpath/libtomar_grouped_residue.dylib", source, "-o", binary], check=True)
        library = ctypes.CDLL(str(binary))
        vec = lambda dtype: np.ctypeslib.ndpointer(dtype=dtype, ndim=1, flags="C_CONTIGUOUS")
        library.fan_probe.argtypes = [ctypes.c_int, vec(np.int64), vec(np.int64), vec(np.int64),
                                     vec(np.int32), vec(np.int32), ctypes.c_int64, ctypes.c_int,
                                     ctypes.c_int, vec(np.int64)]
        goal = 2 * Fraction(TARGET) * SCALE**2
        ceiling = -(-goal.numerator // goal.denominator)
        for file in keys:
            part = int(file.stem[4:])
            record = json.loads(file.read_bytes())
            npz = file.with_suffix(".npz")
            assert record["target"] == TARGET and record["part"] == part
            assert record["remaining_sha256"] == sha(npz)
            dest = args.out / file.name
            with np.load(npz, allow_pickle=False) as data:
                pos = data["pos"]
            assert pos.dtype == np.int32 and pos.ndim == 1 and np.array_equal(pos, np.unique(pos))
            assert len(pos) == record["remaining"]
            codefile = args.certdir / (file.stem + ".npy")
            headerfile = args.certdir / file.name
            assert sha(codefile) == pinned[codefile.name] and sha(headerfile) == pinned[headerfile.name]
            meta = json.loads(headerfile.read_bytes())
            assert meta["hi"] == header["hi"] and meta["lo"] == header["lo"] and meta["part"] == part
            codes = np.load(codefile, allow_pickle=False)
            assert all(0 <= i < len(codes) and int(codes[i]) >= 253 for i in pos)
            rows = []
            if len(pos):
                prefix = np.r_[1, 1 + np.cumsum(np.where(codes < 11, 1, -1))]
                for start in range(0, len(pos), 32):
                    batch = pos[start:start + 32]
                    lo, hi = [], []
                    for i in batch:
                        center, half = leaf_box(prefix, codes, part, int(i), meta)
                        lo.append(scaled([Fraction(x) - Fraction(y) for x, y in zip(center, half)]))
                        hi.append(scaled([Fraction(x) + Fraction(y) for x, y in zip(center, half)]))
                    stats = np.zeros((len(batch), 4), dtype=np.int64)
                    rc = library.fan_probe(len(batch), np.array(lo).ravel(), np.array(hi).ravel(), base,
                                           points.ravel(), counts, ceiling, DEPTH, CALLS, stats.ravel())
                    assert rc == 0
                    rows += [[int(i), int(codes[i]), *map(int, row)] for i, row in zip(batch, stats)]
            result = {"part": part, "target": TARGET, "parent_sha256": sha(file), "positions_sha256": sha(npz),
                      "table_sha256": TABLE, "driver_sha256": sha(Path(__file__)), "source_sha256": sha(source),
                      "decoder_sha256": sha(HERE / "tomar_leaf_box.py"), "library_sha256": sha(binary),
                      "header_sha256": pinned[headerfile.name], "codes_sha256": pinned[codefile.name],
                      "depth": DEPTH, "calls": CALLS, "checked": len(rows),
                      "closed": sum(r[-1] for r in rows), "leaves": rows}
            tmpfile = dest.with_suffix(".partial")
            if dest.exists():
                assert json.loads(dest.read_bytes()) == result
            else:
                tmpfile.write_text(json.dumps(result, separators=(",", ":")) + "\n")
                tmpfile.replace(dest)
            print(part, len(rows), result["closed"], flush=True)

if __name__ == "__main__":
    main()
