"""Certify a dihedral fundamental domain and enclose every possible four-body minimizer."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import numpy as np
from flint import arb, ctx
from replay_four_body_multifan import (
    native,
    check_header,
    core,
    POINT_SHA,
    MATERIALS,
    digest,
)
from make_four_body_material import make, FINE, COORD, rational

CUT = Q(215, 1024)
TARGET = Q(".836490455462")
LIMITS = ((17, 511), (30, 4095))
ROOT_HASH = "526f23db1edc790dc0057503739bb567b57ab47b2e2429959e9899e8dba87fd9"


def stage(a, library, mat, maxdepth, calls):
    pts, hi, counts = mat
    l = np.zeros((len(a), 11), dtype=np.int64)
    r = l.copy()
    l[:, :8], r[:, :8] = a[:, :8], a[:, 8:]
    res = np.zeros((len(a), 4), dtype=np.int64)
    goal = int(2 * TARGET * COORD**2) + 1
    assert (
        library.fan_probe(
            len(a),
            l.ravel(),
            r.ravel(),
            hi,
            pts,
            counts,
            goal,
            maxdepth,
            calls,
            res.ravel(),
        )
        == 0
    )
    ok = res[:, 3].astype(bool)
    assert np.all((res[:, 3] == 0) | (res[:, 3] == 1)) and np.all(
        res[ok, 0] == 2 * res[ok, 2] - 1
    )
    return res, ok


def go(arg):
    tree, work, p = arg
    h = json.loads((tree / "part0000.json").read_bytes())
    boxes, nroots, codes = native(tree, p, h)
    idx = np.flatnonzero(boxes[:, 2] <= int(CUT * FINE))
    active = boxes[idx].copy()
    active[:, 10] = np.minimum(active[:, 10], int(CUT * FINE))
    count = len(idx)
    log = []
    lib, tables = core(work)
    for mat, lim in zip(tables, LIMITS, strict=True):
        res, ok = stage(active, lib, mat, *lim)
        log.append(
            dict(
                inputs=len(active),
                passed=int(ok.sum()),
                checks=int(res[:, 0].sum()),
                input_box_sha=sha256(active.tobytes()).hexdigest(),
                result_sha=sha256(res.tobytes()).hexdigest(),
            )
        )
        active, idx = active[~ok], idx[~ok]
    return dict(
        part=p,
        source_sha=codes,
        original=len(boxes),
        roots=nroots,
        in_fundamental_domain=count,
        stages=log,
        remaining_indices=idx.tolist(),
        remaining_boxes=active.tolist(),
    )


def main():
    ap = ArgumentParser()
    ap.add_argument("tree", type=Path)
    ap.add_argument("work", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--compiler", default="c++")
    a = ap.parse_args()
    tree = a.tree.resolve()
    work = a.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    assert __debug__ and digest(tree / "SHA256SUMS") == ROOT_HASH
    for line in (tree / "SHA256SUMS").read_text().splitlines():
        d, f = line.split()
        assert digest(tree / f) == d
    h = json.loads((tree / "part0000.json").read_bytes())
    ext = check_header(h, TARGET)
    assert arb.pi() / 15 < rational(CUT) and CUT < Q(h["hi"][2])
    here = Path(__file__).parent
    lib = work / "four-body-multifan.so"
    subprocess.run(
        [
            a.compiler,
            "-O3",
            "-std=c++20",
            "-shared",
            "-fPIC",
            str(here / "four_body_multifan.cpp"),
            "-o",
            str(lib),
        ],
        check=True,
    )
    for sha, (n, d) in zip(POINT_SHA, MATERIALS, strict=True):
        tab = make(h, n, d)
        assert sha256(tab["points"].tobytes()).hexdigest() == sha
        np.savez_compressed(work / f"material-{n}.npz", **tab)
    parts = sorted(range(24), key=lambda p: -(tree / f"part{p:04}.npy").stat().st_size)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(go, [(tree, work, p) for p in parts]))
    rows.sort(key=lambda r: r["part"])
    assert (
        sum(r["original"] for r in rows) == 808905
        and sum(r["roots"] for r in rows) == 4**8
    )
    arr = np.asarray([b for r in rows for b in r["remaining_boxes"]], dtype=np.int64)
    report = dict(
        target=str(TARGET),
        cut=str(CUT),
        source_manifest=ROOT_HASH,
        source_sha={
            k: digest(here / k)
            for k in (
                "localize_four_body_global.py",
                "replay_four_body_multifan.py",
                "make_four_body_material.py",
                "four_body_multifan.cpp",
            )
        },
        actual_library_sha=digest(lib),
        compiler=subprocess.check_output(
            [a.compiler, "--version"], text=True
        ).splitlines()[0],
        table_hashes=POINT_SHA,
        denominator=FINE,
        exterior=ext,
        original=808905,
        active_roots=sum(r["in_fundamental_domain"] for r in rows),
        remaining=len(arr),
        coordinate_hull=arr[:, :8].min(0).tolist() + arr[:, 8:].max(0).tolist(),
        parts=rows,
    )
    a.out.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    print(report["active_roots"], report["remaining"], report["coordinate_hull"])


if __name__ == "__main__":
    main()
