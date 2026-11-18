"""Independently prove a four-body lower bound using only the original split topology."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import islice, product
from pathlib import Path
import ctypes
import json
import os
import platform
import shutil
import flint
import subprocess
import sys
import numpy as np
from flint import arb, ctx
from make_four_body_material import COORD, FINE, make, rational

if not __debug__:
    raise SystemExit("Run with assertions enabled.")
HERE = Path(__file__).parent
MATERIALS = ((24, 96), (64, 384))
POINT_SHA = (
    "119f6357493a034dca1ce0dc4be5229ff2b625701403d0065725eb90d1c00a99",
    "0ad9cbb26f37e730732d7f41c4c61a36ba49516a637318a721b27159f94fa7f3",
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def scaled_endpoints(v):
    exact = [F(x) * FINE for x in v]
    assert all(y.denominator == 1 for y in exact)
    return np.array([int(y) for y in exact], dtype=np.int64)


def check_header(h, target):
    ctx.prec = 128
    assert h["ns"] == [3, 5, 7] and h["split0"] == [4] * 8 and h["nparts"] == 24
    assert len(h["lo"]) == len(h["hi"]) == 8
    lo, hi = scaled_endpoints(h["lo"]), scaled_endpoints(h["hi"])
    assert np.all(lo < hi) and np.all(np.abs(lo) < FINE) and np.all(np.abs(hi) < FINE)
    for j, n in ((2, 5), (5, 7)):
        assert (
            lo[j] == 0 and rational(h["hi"][j]) >= arb.pi() * (1 if j == 2 else 2) / n
        )
    exterior = []
    for n, j in ((3, 0), (5, 3), (7, 6)):
        r = 1 - 1 / (2 * (arb.pi() / (2 * n)).cos())
        w = min(-F(h["lo"][j]), -F(h["lo"][j + 1]), F(h["hi"][j]), F(h["hi"][j + 1]))
        bound = arb.pi() * (arb(1) / 4 + r * r) / 2 + rational(w) * (arb(1) / 2 + r)
        assert w > 0 and r > 0 and bound > rational(target)
        exterior.append(str(bound))
    return exterior


def native(tree, p, header):
    tag = f"part{p:04d}"
    meta = json.loads((tree / (tag + ".json")).read_bytes())
    assert meta["part"] == p and all(
        meta[k] == header[k] for k in ("ns", "lo", "hi", "split0", "nparts")
    )
    path = tree / (tag + ".npy")
    code = np.load(path, allow_pickle=False)
    assert code.ndim == 1 and code.dtype == np.uint8
    l = scaled_endpoints(header["lo"])
    b = scaled_endpoints(header["hi"])
    cw = (b - l) // 4
    assert (4 * cw == b - l).all()
    a = np.empty((np.count_nonzero(code == 254), 16), dtype=np.int64)
    pos = leaf = roots = 0
    for cell in islice(product(range(4), repeat=8), p, None, 24):
        roots += 1
        lo = l + cw * cell
        stack = [(lo, lo + cw)]
        while stack:
            lo, hi = stack.pop()
            axis = int(code[pos])
            pos += 1
            if axis == 254:
                a[leaf, :8] = lo
                a[leaf, 8:] = hi
                leaf += 1
                continue
            assert 0 <= axis < 8 and (lo[axis] + hi[axis]) % 2 == 0
            mid = (lo[axis] + hi[axis]) // 2
            sublo = lo.copy()
            sublo[axis] = mid
            subhi = hi.copy()
            subhi[axis] = mid
            stack.extend(((sublo, hi), (lo, subhi)))
    assert pos == len(code) and leaf == len(a)
    return a, roots, digest(path)


@lru_cache(None)
def core(work):
    dll = ctypes.CDLL(str(work / "four-body-multifan.so"))
    vec = lambda t: np.ctypeslib.ndpointer(dtype=t, ndim=1, flags="C_CONTIGUOUS")
    dll.fan_probe.argtypes = [
        ctypes.c_int,
        vec(np.int64),
        vec(np.int64),
        vec(np.int64),
        vec(np.int32),
        vec(np.int32),
        ctypes.c_int64,
        ctypes.c_int,
        ctypes.c_int,
        vec(np.int64),
    ]
    params = []
    for j, (sub, disk) in enumerate(MATERIALS):
        with np.load(work / f"material-{sub}.npz", allow_pickle=False) as src:
            pts = src["points"]
            counts = src["counts"]
            hi = src["original_hi"]
        assert (
            counts.tolist() == [disk, 3 * sub, 5 * sub, 7 * sub, 0]
            and counts.dtype == np.int32
        )
        assert pts.dtype == np.int32 and pts.shape == (1025, sum(counts), 2)
        assert sha256(pts.tobytes()).hexdigest() == POINT_SHA[j]
        params.append((np.ascontiguousarray(pts).ravel(), hi, counts))
    return dll, params


def part(job):
    tree, work, p, target, limits = job
    header = json.loads((tree / "part0000.json").read_bytes())
    boxes, roots, sha = native(tree, p, header)
    original = np.arange(len(boxes), dtype=np.int64)
    total = len(boxes)
    dll, params = core(work)
    exact = 2 * target * COORD**2
    # Nonstrict integer comparison is sufficient for a strict real threshold.
    goal = exact.numerator // exact.denominator + 1
    stages = []
    for (depth, calls), (pts, hi, counts) in zip(limits, params, strict=True):
        q = boxes[original]
        l = np.zeros((len(q), 11), dtype=np.int64)
        r = l.copy()
        l[:, :8] = q[:, :8]
        r[:, :8] = q[:, 8:]
        assert np.array_equal(hi[:8], scaled_endpoints(header["hi"])) and hi[8] > 0
        result = np.zeros((len(q), 4), dtype=np.int64)
        assert (
            dll.fan_probe(
                len(q),
                l.ravel(),
                r.ravel(),
                np.ascontiguousarray(hi),
                pts,
                counts,
                goal,
                depth,
                calls,
                result.ravel(),
            )
            == 0
        )
        passed = result[:, 3] == 1
        assert np.all(passed | (result[:, 3] == 0))
        assert np.all(result[passed, 0] == 2 * result[passed, 2] - 1)
        input_ids_sha = sha256(original.tobytes()).hexdigest()
        original = original[~passed]
        stages.append(
            dict(
                inputs=len(q),
                input_indices_sha256=input_ids_sha,
                input_box_endpoints_sha256=sha256(q.tobytes()).hexdigest(),
                integer_result_ledger_sha256=sha256(result.tobytes()).hexdigest(),
                remaining_indices_sha256=sha256(original.tobytes()).hexdigest(),
                passed=int(passed.sum()),
                remaining=len(original),
                checks=int(result[:, 0].sum()),
                max_depth=int(result[:, 1].max(initial=0)),
            )
        )
    return dict(
        part=p,
        root_cells=roots,
        source_sha256=sha,
        source_meta_sha256=digest(tree / f"part{p:04d}.json"),
        native_box_endpoints_sha256=sha256(boxes.tobytes()).hexdigest(),
        native_leaves=total,
        stages=stages,
        remaining_native_indices=original.tolist(),
    )


def main():
    a = ArgumentParser()
    a.add_argument("tree", type=Path)
    a.add_argument("work", type=Path)
    a.add_argument("output", type=Path)
    a.add_argument("--target", default=".8364")
    a.add_argument("--workers", type=int, default=4)
    a.add_argument("--depths", type=int, nargs=2, default=[23, 25])
    a.add_argument("--calls", type=int, nargs=2, default=[4095, 4095])
    a.add_argument("--compiler", default="c++")
    opt = a.parse_args()
    tree, work = opt.tree.resolve(), opt.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    target = F(opt.target)
    h = json.loads((tree / "part0000.json").read_bytes())
    manifest = tree / "SHA256SUMS"
    assert (
        digest(manifest)
        == "526f23db1edc790dc0057503739bb567b57ab47b2e2429959e9899e8dba87fd9"
    )
    entries = {
        p: sha for sha, p in (s.split() for s in manifest.read_text().splitlines())
    }
    assert set(entries) == {
        f"part{i:04d}.{ext}" for i in range(24) for ext in ("json", "npy")
    }
    assert all(digest(tree / name) == sha for name, sha in entries.items())
    exterior = check_header(h, target)
    assert (
        opt.workers >= 1
        and all(0 <= d <= 48 for d in opt.depths)
        and all(c >= 1 for c in opt.calls)
    )
    src = HERE / "four_body_multifan.cpp"
    env = dict(
        os.environ,
        CLANG_MODULE_CACHE_PATH=os.environ.get(
            "CLANG_MODULE_CACHE_PATH", "/Volumes/data/CodexTemp/cache/clang-modules"
        ),
    )
    subprocess.run(
        [
            opt.compiler,
            "-O3",
            "-std=c++20",
            "-shared",
            "-fPIC",
            str(src),
            "-o",
            str(work / "four-body-multifan.so"),
        ],
        check=True,
        env=env,
    )
    for j, (sub, disk) in enumerate(MATERIALS):
        data = make(h, sub, disk)
        assert sha256(data["points"].tobytes()).hexdigest() == POINT_SHA[j]
        np.savez_compressed(work / f"material-{sub}.npz", **data)
    limits = list(zip(opt.depths, opt.calls))
    parts = sorted(range(24), key=lambda p: -(tree / f"part{p:04d}.npy").stat().st_size)
    from concurrent.futures import as_completed

    with ProcessPoolExecutor(max_workers=opt.workers) as pool:
        jobs = [pool.submit(part, (tree, work, p, target, limits)) for p in parts]
        rows = []
        for result in as_completed(jobs):
            row = result.result()
            rows.append(row)
            print(
                row["part"],
                row["native_leaves"],
                len(row["remaining_native_indices"]),
                flush=True,
                file=sys.stderr,
            )
    rows.sort(key=lambda r: r["part"])
    assert len(rows) == 24 and sum(r["root_cells"] for r in rows) == 4**8
    remain = sum(len(r["remaining_native_indices"]) for r in rows)
    totals = [
        {
            k: sum(r["stages"][t][k] for r in rows)
            for k in ("inputs", "passed", "remaining", "checks")
        }
        for t in range(2)
    ]
    final = dict(
        target=str(target),
        upstream_tree_manifest_sha256=digest(manifest),
        source_sha256={
            f: digest(HERE / f)
            for f in (
                "replay_four_body_multifan.py",
                "four_body_multifan.cpp",
                "make_four_body_material.py",
            )
        },
        runtime=dict(
            python=platform.python_version(),
            python_flint=flint.__version__,
            numpy=np.__version__,
            compiler=shutil.which(opt.compiler),
            compiler_version=subprocess.check_output(
                [opt.compiler, "--version"], text=True
            ).splitlines()[0],
            compiler_flags=["-O3", "-std=c++20", "-shared", "-fPIC"],
            compiled_library_sha256=digest(work / "four-body-multifan.so"),
        ),
        material_archive_sha256=[
            digest(work / f"material-{s}.npz") for s, d in MATERIALS
        ],
        complete=(remain == 0),
        native_leaves=sum(r["native_leaves"] for r in rows),
        initial_cells=4**8,
        material_point_hashes=list(POINT_SHA),
        limits=limits,
        exterior=exterior,
        stages=totals,
        remaining=remain,
        parts=rows,
    )
    opt.output.write_text(json.dumps(final, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: final[k]
                for k in ("target", "complete", "native_leaves", "remaining", "stages")
            }
        )
    )
    if remain:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
