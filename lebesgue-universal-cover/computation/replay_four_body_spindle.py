"""Replay a whole complement of a dyadic window from all 287 normalized boxes."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import platform
import numpy as np
import scipy
import flint
from four_body_spindle_solver import P, S, vertices, work

HERE = Path(__file__).resolve().parent
MATERIAL = HERE / "four-body-spindle-material.npz"
BASE = HERE / "four-body-normalized-localization.json"
SOURCE_SHA = "7293b6f161d1b2699d396d43457e6a6e10ee33f0428db286f6c11f4cdf1d0051"
RECEIPT_SHA = "5f76843fe5c961f11e0e1ae2dddfd95fe40c4e5666110c36e4cfdc677e993721"
# Each endpoint is in units of 2^-44.
LO = [
    -46648066048,
    -35770466304,
    486822051840,
    439573938176,
    -59545944064,
    12324394172416,
    49037967360,
    -225511604224,
]
HI = [
    -43129569280,
    -32252035072,
    522006429696,
    443092369408,
    -56027512832,
    12359578550272,
    52556398592,
    -221993172992,
]
ORDER = [2, 5, 0, 1, 3, 4, 6, 7]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify_material():
    data = np.load(MATERIAL)
    counts = data["counts"]
    pts = data["points"].astype(object)
    assert counts.tolist() == [205, 152, 210, 253, 0] and pts.shape == (sum(counts), 2)
    assert np.all(np.sum(np.abs(pts), axis=1) <= P)
    off = 0
    for b, n in enumerate([0, 3, 5, 7]):
        arr = pts[off : off + counts[b]]
        off += counts[b]
        if n:
            vv = vertices(n, 0).astype(object)
            # Vertex approximation error is < sqrt(2)/P.
            diff = arr[:, None, :] - vv
            assert np.all(np.sum(diff * diff, axis=2) < (P - 7) ** 2)
        else:
            assert np.all(np.sum(arr * arr, axis=1) < (P // 2 - 5) ** 2)
    return digest(MATERIAL)


def exterior(rows, low, high):
    roots = []
    for idx, box in rows:
        a = box[:8].copy()
        b = box[8:].copy()
        for k in ORDER:
            if a[k] < low[k]:
                c = b.copy()
                c[k] = min(b[k], low[k])
                roots.append((idx, -1 - k, a + c))
                a[k] = c[k]
            if b[k] > high[k]:
                c = a.copy()
                c[k] = max(a[k], high[k])
                roots.append((idx, 1 + k, c + b))
                b[k] = c[k]
            if a[k] == b[k]:
                break
    return roots


def main():
    ap = ArgumentParser()
    ap.add_argument(
        "geometry", type=Path, help="v9 independent/gen8.py (proposal generation only)"
    )
    ap.add_argument("workdir", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--compiler", default="c++")
    a = ap.parse_args()
    assert (
        __debug__ and digest(BASE) == RECEIPT_SHA and digest(a.geometry) == SOURCE_SHA
    )
    a.workdir = a.workdir.resolve()
    a.workdir.mkdir(parents=True, exist_ok=True)
    os.environ["FOUR_SPINDLE_GEOMETRY"] = str(a.geometry.resolve().parent)
    cpp = a.workdir / "four_body_spindle.cpp"
    with (HERE / "four_body_spindle.patch").open("rb") as inp:
        subprocess.run(
            ["patch", "-s", "-o", str(cpp), str(HERE / "four_body_multifan.cpp")],
            stdin=inp,
            check=True,
        )
    lib = a.workdir / "four_body_spindle.so"
    flags = [
        "-O3",
        "-std=c++20",
        "-shared",
        "-fPIC",
        "-fsanitize=signed-integer-overflow",
        "-fsanitize-undefined-trap-on-error",
    ]
    subprocess.run([a.compiler, *flags, str(cpp), "-o", str(lib)], check=True)
    os.environ["FOUR_SPINDLE_LIBRARY"] = str(lib)
    src = json.loads(BASE.read_text())
    assert src["denominator"] == S and src["remaining"] == 287
    p = src["parts"][17]
    assert all(v["part"] == 17 or not v["remaining_boxes"] for v in src["parts"])
    rows = list(zip(p["remaining_indices"], p["remaining_boxes"], strict=True))
    assert len(rows) == 287
    stages = []
    report = dict(
        source_receipt=RECEIPT_SHA,
        geometry=SOURCE_SHA,
        table=verify_material(),
        denominator=S,
        low=LO,
        high=HI,
        compiled_source_sha=digest(cpp),
        actual_library_sha=digest(lib),
        flags=flags,
        compiler=subprocess.check_output(
            [a.compiler, "--version"], text=True
        ).splitlines()[0],
        versions=[
            platform.python_version(),
            np.__version__,
            scipy.__version__,
            flint.__version__,
        ],
        sources={
            name: digest(HERE / name)
            for name in (
                "four_body_multifan.cpp",
                "four_body_spindle.patch",
                "four_body_spindle_solver.py",
                "replay_four_body_spindle.py",
            )
        },
        stages=stages,
    )
    jobs = exterior(rows, LO, HI)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for budget, depth in ((255, 40), (511, 48), (1023, 48), (2047, 64)):
            if not jobs:
                break
            pending = []
            total = 0
            hashes = sha256()
            for arg, r in zip(
                jobs, ex.map(work, [(*v, budget, depth) for v in jobs]), strict=True
            ):
                assert r["id"] == arg[0] and r["label"] == arg[1]
                hashes.update(
                    (
                        json.dumps([arg, r], sort_keys=True, separators=(",", ":"))
                        + "\n"
                    ).encode()
                )
                pending.extend((r["id"], r["label"], v) for v in r["remaining"])
                total += r["vis"]
            stages.append(
                dict(
                    inputs=len(jobs),
                    budget=budget,
                    depth=depth,
                    visits=total,
                    left=len(pending),
                    stage_sha=hashes.hexdigest(),
                )
            )
            print(stages[-1], flush=True)
            jobs = pending
    report["remaining"] = jobs
    a.output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    assert not jobs, "Some of the complement remains"


if __name__ == "__main__":
    main()
