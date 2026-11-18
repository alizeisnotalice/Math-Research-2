"""Recover the certified basin complement from the globally localized spindle window."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from replay_four_body_spindle import HERE, SOURCE_SHA, digest, verify_material
from four_body_spindle_solver import S, work

WINDOW_SHA = "4ff0cbf100a36280ec5fe5d717ef639b3f11e658f5243104eb18aea9e402d171"
LOW = [
    -45240680448,
    -34363080704,
    502655025152,
    440981323776,
    -58138558464,
    12340227145728,
    50445287424,
    -224104218624,
]
HIGH = [
    -44536954880,
    -33659355136,
    506173456384,
    441684983808,
    -57434898432,
    12343745576960,
    51149012992,
    -223400558592,
]


def seed_parts(lower, upper):
    a = np.asarray(lower, dtype=np.int64).copy()
    b = np.asarray(upper, dtype=np.int64).copy()
    scales = np.array([1, 1, 0.16, 1, 1, 0.11, 1, 1])
    slabs = []
    for k in [0, 1, 2, 5, 3, 4, 6, 7]:
        assert a[k] < LOW[k] < HIGH[k] < b[k]
        c = b.copy()
        c[k] = LOW[k]
        slabs.append((-k - 1, np.r_[a, c]))
        a[k] = LOW[k]
        c = a.copy()
        c[k] = HIGH[k]
        slabs.append((1 + k, np.r_[c, b]))
        b[k] = HIGH[k]
    out = []
    for label, box in slabs:
        todo = [box]
        while len(todo) < 8:
            spans = np.array([(r[8:] - r[:8]) * scales for r in todo])
            i, k = np.unravel_index(np.argmax(spans), spans.shape)
            r = todo.pop(i)
            lo = r.copy()
            hi = r.copy()
            mid = (r[k] + r[8 + k]) // 2
            lo[k + 8] = mid
            hi[k] = mid
            assert r[k] < mid < r[k + 8]
            todo.extend([lo, hi])
        out.extend((label, r.tolist()) for r in todo)
    return [(j, lab, box) for j, (lab, box) in enumerate(out)]


def main():
    ap = ArgumentParser()
    ap.add_argument("geometry", type=Path)
    ap.add_argument("workdir", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--compiler", default="c++")
    args = ap.parse_args()
    assert __debug__ and digest(args.geometry) == SOURCE_SHA
    src = HERE / "four-body-spindle-window.json"
    assert digest(src) == WINDOW_SHA
    prior = json.loads(src.read_text())
    assert prior["denominator"] == S and not prior["remaining"]
    for name, sha in prior["sources"].items():
        assert digest(HERE / name) == sha
    assert verify_material() == prior["table"]
    args.workdir = args.workdir.resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)
    cpp = args.workdir / "four_body_spindle.cpp"
    lib = args.workdir / "four_body_spindle.so"
    with (HERE / "four_body_spindle.patch").open("rb") as patch:
        subprocess.run(
            ["patch", "-s", "-o", str(cpp), str(HERE / "four_body_multifan.cpp")],
            stdin=patch,
            check=True,
        )
    assert digest(cpp) == prior["compiled_source_sha"]
    subprocess.run(
        [args.compiler, *prior["flags"], str(cpp), "-o", str(lib)], check=True
    )
    os.environ["FOUR_SPINDLE_LIBRARY"] = str(lib)
    os.environ["FOUR_SPINDLE_GEOMETRY"] = str(args.geometry.resolve().parent)
    jobs = seed_parts(prior["low"], prior["high"])
    stages = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for budget, depth in [
            (511, 48),
            (1023, 48),
            (2047, 48),
            (4095, 64),
            (8191, 64),
            (16383, 64),
        ]:
            if not jobs:
                break
            h = sha256()
            pending = []
            total = 0
            for arg, r in zip(
                jobs,
                ex.map(work, [(*item, budget, depth) for item in jobs]),
                strict=True,
            ):
                assert r["id"] == arg[0] and r["label"] == arg[1]
                h.update(
                    (
                        json.dumps([arg, r], sort_keys=True, separators=(",", ":"))
                        + "\n"
                    ).encode()
                )
                total += r["vis"]
                pending.extend((r["id"], r["label"], b) for b in r["remaining"])
            stages.append(
                dict(
                    inputs=len(jobs),
                    budget=budget,
                    depth=depth,
                    visits=total,
                    left=len(pending),
                    sha=h.hexdigest(),
                )
            )
            print(stages[-1], flush=True)
            jobs = pending
    report = dict(
        prior=WINDOW_SHA,
        source=digest(Path(__file__)),
        target="0.836490455462",
        denominator=S,
        low=LOW,
        high=HIGH,
        stages=stages,
        remaining=jobs,
        compiler=subprocess.check_output(
            [args.compiler, "--version"], text=True
        ).splitlines()[0],
        compiled_source_sha=digest(cpp),
        actual_library_sha=digest(lib),
        flags=prior["flags"],
    )
    args.output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    assert not jobs, "Some of the basin complement remains"


if __name__ == "__main__":
    main()
