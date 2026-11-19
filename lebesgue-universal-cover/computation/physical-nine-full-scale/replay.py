"""Replay the physical nine-input lower bound and the fixed-placement upper tests."""
from argparse import ArgumentParser
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
PACKAGES = {
    "scale.zip": "3c22553b446386cb416c10d7505fc93e449feac4ad84defa075bab0368131542",
    "poses.zip": "488ab6e7e6f60de999c676fb1e833c135ec800ff79ef468bf78cc70f40f50d72",
    "additions.zip": "92b068d9b5e9ed2f4007506404fa7c2ba217fe5804470b25557f054d6c11c485",
    "point-counts.zip": "3f65d7629c9dddb03fc875241186d552174f419ca5721aab5409995556930827",
    "curve-cross.zip": "f39bc5d859565c6a907f90ae5ec0ffd008a977472b4573748695d1935a402384",
    "curve-points.zip": "f43d11c27248440fb7c41871605f3dbb071c4758c9d3e09473f5f88ea22dcace",
}


def read(file):
    return json.loads(file.read_bytes())


def run(folder, script, *args):
    subprocess.run([sys.executable, "-I", "-B", script, *map(str, args)], cwd=folder, check=True)


def equal(expected, actual):
    assert expected.read_bytes() == actual.read_bytes(), expected.name


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("clang_cache", type=Path)
    ap.add_argument("--r286-receipt", type=Path, help="receipt.json from the R286 basins-replay directory")
    args = ap.parse_args()
    assert __debug__
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    for name, expected in PACKAGES.items():
        path = HERE / name
        assert sha256(path.read_bytes()).hexdigest() == expected
        with ZipFile(path) as z:
            assert all(not Path(x).is_absolute() and ".." not in Path(x).parts for x in z.namelist())
            z.extractall(work / path.stem)

    scaled = work / "scaled"
    run(work / "scale", "replay.py", HERE.parent, scaled)
    proof = read(scaled / "final/result.json")
    G = Fraction(120697911165410052, 2**57)
    scale = Fraction(proof["sigma"])
    assert scale == Fraction(proof["s"]) == Fraction(499999503, 500000000)
    assert G / scale**2 == Fraction(proof["area"])
    assert sum(proof["tables"][0]["tested"]) == 119208
    points = {}
    for n in [6131, 9769]:
        d = scaled / f"finite{n}"
        assert gzip.decompress((HERE / f"finite{n}.json.gz").read_bytes()) == gzip.decompress((d / "points.json.gz").read_bytes())
        row = read(d / "receipt.json")
        assert Fraction(row["general_lower"]) == G / scale**2 and row["N"] == n
        points[str(n)] = {"decimal": row["decimal"], "columns": row["columns"], "points_sha": row["points_sha"]}
    hulls = work / "hulls.json"
    limits = work / "limits.json"
    src = work / "point-counts"
    run(src, "extremes51.py", scaled / "final/physical.json", hulls)
    run(src, "finite_limits.py", limits)
    equal(src / "hulls51.json", hulls)
    equal(src / "finite-limits.json", limits)
    sizes = [len(i["points"]) for i in read(hulls)]
    assert sizes == [114, 11, 11, 15]

    poses = work / "poses"
    specific = {}
    for key, initial, reference, n in [
            ("nine", "refine-search.json", "best-nine-v2-exact.json", 131072),
            ("six", "ref-six-search.json", "reference-six-fine.json", 262144),
            ("seven", "ref-seven-search.json", "ref-seven-exact.json", 32768)]:
        out = work / (key + "-pose.json")
        poly = work / (key + "-host.json")
        run(poses, "verify_upper.py", poses / "inputs", poses / "legacy", poses / initial,
            poly, out, "--samples", n, *(["--diagnostics"] if key == "nine" else []))
        equal(poses / reference, out)
        d = read(out)
        specific[key] = {"area_upper": d["exact_outer_area"]}
        if key == "nine":
            specific[key]["area_lower"] = d["proved_inner_area"]
    equal(poses / "best-nine-v2-host.json", work / "nine-host.json")
    extra = work / "additions"
    run(extra, "review_extension.py", poses / "inputs/p128.json", work / "nine-host.json",
        work / "nine-pose.json", extra / "p128-extension-facets.json", work / "extra.json")
    equal(extra / "p128-fits-audit.json", work / "extra.json")

    cross = work / "curve-cross/nine-w-r7-cross"
    run(cross, "verify.py", work / "cross", args.clang_cache.resolve())
    whole = read(work / "cross/receipt.json")
    assert whole["v2"]["host_sha256"] == read(work / "nine-pose.json")["host_sha256"]
    scope = []
    if args.r286_receipt:
        scope = ["--r286-actual-replay-receipt", args.r286_receipt.resolve()]
    run(work / "curve-points/nine-w-finite", "verify.py", cross, work / "cross/receipt.json",
        work / "finite-W.json", *scope)
    finite_w = read(work / "finite-W.json")
    assert finite_w["status"] == "PASS"
    W = whole["v1"]["W"]
    result = {"unrestricted_lower": proof["area"], "decimal": proof["decimal"], "scale": proof["sigma"],
              "new_P_angular_instances": proof["tables"][0]["tested"],
              "angular_row_surplus": proof["tables"][0]["decimal_projection_spare"],
              "physical_hull_sizes": sizes, "finite_lower_bounds": points,
              "point_method_upper_bounds": read(limits), "specific_placements": specific,
              "W_added_area_v1": [W["added_area_lower"], W["added_area_upper"]],
              "two_R7_v2_joint_area_increment_upper": whole["v2"]["simultaneous_true_area_increment_upper"],
              "W_cloud_points": [len(i["points"]) for i in finite_w["clouds"]],
              "fixed_six_checked": bool(args.r286_receipt)}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
