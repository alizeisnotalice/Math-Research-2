"""Close or report all remaining four-body leaves using the frozen C++ kernel with denser material."""

from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import numpy as np
from replay_four_body_multifan import HERE, check_header, digest, native
from make_four_body_material import COORD, make

EXPECTED = "8227769fe2abf634a43ef2975cf45366d2802eef40d781dd9b1674cc15a6af4f"


def finish(tree, previous, work, output, compiler, depth, calls):
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    work.mkdir(parents=True, exist_ok=True)
    previous_data = json.loads(previous.read_bytes())
    target = F(previous_data["target"])
    h = json.loads((tree / "part0000.json").read_bytes())
    check_header(h, target)
    assert digest(tree / "SHA256SUMS") == previous_data["upstream_tree_manifest_sha256"]
    for name, code_sha in previous_data["source_sha256"].items():
        assert digest(HERE / name) == code_sha
    assert (
        previous_data["initial_cells"] == 4**8
        and previous_data["native_leaves"] == 808905
    )
    assert [r["part"] for r in previous_data["parts"]] == list(range(24))
    assert 0 <= depth <= 48 and 1 <= calls < 2**24
    data = make(h, 128, 512)
    pt = data["points"]
    assert sha256(pt.tobytes()).hexdigest() == EXPECTED
    vec = lambda t: np.ctypeslib.ndpointer(dtype=t, ndim=1, flags="C_CONTIGUOUS")
    library = work / "four-body-finish.so"
    flags = ["-O3", "-std=c++20", "-shared", "-fPIC"]
    subprocess.run(
        [compiler, *flags, str(HERE / "four_body_multifan.cpp"), "-o", str(library)],
        check=True,
    )
    dll = ctypes.CDLL(str(library))
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
    threshold = 2 * target * COORD**2
    goal = threshold.numerator // threshold.denominator + 1
    rows = []
    for prev in previous_data["parts"]:
        p = prev["part"]
        ix = np.asarray(prev["remaining_native_indices"], dtype=np.int64)
        assert ix.tolist() == sorted(set(ix.tolist())) and (ix >= 0).all()
        assert len(ix) == prev["stages"][-1]["remaining"]
        assert (
            sha256(ix.tobytes()).hexdigest()
            == prev["stages"][-1]["remaining_indices_sha256"]
        )
        if not len(ix):
            continue
        b, roots, code_sha = native(tree, p, h)
        assert prev["root_cells"] == roots and prev["source_sha256"] == code_sha
        assert sha256(b.tobytes()).hexdigest() == prev["native_box_endpoints_sha256"]
        a = b[ix]
        lo = np.zeros((len(a), 11), dtype=np.int64)
        hi = lo.copy()
        lo[:, :8] = a[:, :8]
        hi[:, :8] = a[:, 8:]
        result = np.zeros((len(a), 4), dtype=np.int64)
        assert not dll.fan_probe(
            len(a),
            lo.ravel(),
            hi.ravel(),
            data["original_hi"],
            pt.ravel(),
            data["counts"],
            goal,
            depth,
            calls,
            result.ravel(),
        )
        passed = result[:, 3] == 1
        assert np.all(passed | (result[:, 3] == 0)) and np.all(
            result[passed, 0] == 2 * result[passed, 2] - 1
        )
        rows.append(
            dict(
                part=p,
                source_sha256=code_sha,
                inputs=ix.tolist(),
                input_boxes_sha256=sha256(a.tobytes()).hexdigest(),
                result_ledger=result.tolist(),
                remaining=ix[~passed].tolist(),
            )
        )
    rem = sum(len(r["remaining"]) for r in rows)
    inp = sum(len(r["inputs"]) for r in rows)
    assert inp == previous_data["remaining"]
    report = dict(
        target=str(target),
        previous_sha256=digest(previous),
        material_point_sha256=EXPECTED,
        verifier_sha256=digest(Path(__file__)),
        cpp_sha256=digest(HERE / "four_body_multifan.cpp"),
        compiler_version=subprocess.check_output(
            [compiler, "--version"], text=True
        ).splitlines()[0],
        flags=flags,
        actual_library_sha256=digest(library),
        depth=depth,
        calls=calls,
        all_original=previous_data["native_leaves"],
        inputs=inp,
        remaining=rem,
        checks=sum(v[0] for r in rows for v in r["result_ledger"]),
        complete=(rem == 0),
        parts=rows,
    )
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        {k: report[k] for k in ("target", "inputs", "remaining", "checks", "complete")}
    )
    if rem:
        raise SystemExit(1)


if __name__ == "__main__":
    ap = ArgumentParser()
    for name in ("tree", "previous", "work", "output"):
        ap.add_argument(name, type=Path)
    ap.add_argument("--compiler", default="c++")
    ap.add_argument("--depth", type=int, default=40)
    ap.add_argument("--calls", type=int, default=65535)
    a = ap.parse_args()
    finish(
        a.tree.resolve(),
        a.previous,
        a.work.resolve(),
        a.output,
        a.compiler,
        a.depth,
        a.calls,
    )
