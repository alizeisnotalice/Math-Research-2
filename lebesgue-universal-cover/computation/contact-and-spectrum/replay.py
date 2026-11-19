"""Reproduce the finite moment and rotation examples in the contact compression theorem."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
PACKAGES = {
    "generic": "b95ca0e8123b1b8ad845c4accff4753140196b3b0efcfcb273e6c8d905a52132",
    "explicit": "1ffd265c68b2571962eb0a28ac35ce3a2137aa4d2c9c9c84f06f0c491c04dc6d",
    "rotations": "dee23223b5d2fba485a06c1160807da4dc473f9178ce6f400bfa735f1854b46e",
}


def run(folder, name, *arguments, python=sys.executable):
    subprocess.run([str(python), "-I", "-B", name, *map(str, arguments)], cwd=folder, check=True)


def reproduce(folder, source, destination, expected, *args, parameters_first=False):
    params = [*args, destination] if parameters_first else [destination, *args]
    run(folder, source, *params)
    assert destination.read_bytes() == (folder / expected).read_bytes(), expected
    return json.loads(destination.read_bytes())


def error(n):
    c = (arb.pi() / (2 * n)).cos()
    return 2 * (1 - c) / ((3 + c*c).sqrt() + 1 + c)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--symbolic-python", type=Path, help="Python with SymPy for independent exact integration")
    args = parser.parse_args()
    assert __debug__
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    for name, digest in PACKAGES.items():
        arc = HERE / (name + ".zip")
        assert hashlib.sha256(arc.read_bytes()).hexdigest() == digest
        with ZipFile(arc) as z:
            assert all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in z.namelist())
            z.extractall(work / name)
    g = work / "generic"
    e = work / "explicit"
    o = work / "rotations"
    reproduce(g, "moment_rank.py", work / "rank.json", "moment-rank.json")
    reproduce(g, "circle_examples.py", work / "disks.json", "circle-examples.json")
    reproduce(e, "construct.py", work / "six-windows.json", "receipt.json")
    larger = reproduce(e, "general.py", work / "n3m5.json", "general-3-5.json", 3, 5, parameters_first=True)
    run(e, "boundary_audit.py")
    if args.symbolic_python:
        run(e, "sympy_audit.py", python=args.symbolic_python.absolute())
        run(HERE, "symbolic.py", python=args.symbolic_python.absolute())
    rows = [reproduce(o, "octagon_analytic.py", work / f"angles-{n}.json", expected, n)
            for n, expected in [(8, "octagon_analytic_8.json"), (16, "octagon_analytic.json")]]
    ctx.prec = 192
    A = arb(".84401301378831905312339351556880")
    B = arb(".844082046521218581802948589770767704354")
    bounds = {}
    for n in [173, 174]:
        eps = error(n)
        upper = A + 4*A*eps + arb.pi()*eps**2
        assert (upper > B) if n == 173 else (upper < B)
        bounds[str(n)] = str(upper)
    thresholds = {}
    for n, gap in [(249, ".00001"), (2484, ".0000001")]:
        assert error(n) < arb(gap) < error(n-1)
        thresholds[gap] = {"grid": n, "arcs_through_five": 2*n+5}
    assert larger["d"] == 8 and larger["full_rank"] == 10
    result = {"six_window_certificate": hashlib.sha256((work / "six-windows.json").read_bytes()).hexdigest(),
              "eight_dimensional_rank": larger["full_rank"], "degrees_of_freedom": larger["d"],
              "octagon_whole_angle_margins": [x["interval_minimum"] for x in rows],
              "steiner_area_upper_formula": bounds, "full_pose_gap_thresholds": thresholds}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
