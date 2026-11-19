"""Rebuild the two-cap area obstruction, the pentagon test and the repaired polygons."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import zipfile

from flint import arb, ctx

HERE = Path(__file__).resolve().parent


def unpack(path, target, expected):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    with zipfile.ZipFile(path) as container:
        assert all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in container.namelist())
        container.extractall(target)


def run(source, destination, *options):
    subprocess.run([sys.executable, "-I", "-B", str(source / "replay.py"),
                    str(destination), *options], cwd=source, check=True)


def read(directory, name):
    return json.loads((directory / name).read_bytes())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("--rebuild-base", action="store_true",
                    help="rebuild the C++17 pentagon sweep and all 416676 closed rectangles")
    args = ap.parse_args()
    assert __debug__
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    weighted = work / "weight"
    neighborhood = work / "tube"
    unpack(HERE / "weighted.zip", weighted, "002861922b6be9794c5fb7bd1870305af6e15297f438671793ea3bf2357bc9c9")
    unpack(HERE / "neighborhoods.zip", neighborhood, "1e2754b09ae53b4a6dc5d2a4986133e64c3cd2889ea73126bec5f3a0dfdeadb4")
    run(weighted, work / "obstruction")
    run(neighborhood, work / "neighborhood")
    joint = work / "joint"
    unpack(work / "neighborhood/previous/prior.zip", joint, "9ba5a26d6de13402b28bead83a7627c8f700ae400271b2a543dcde757f50ba78")
    run(joint, work / "first", *(["--base"] if args.rebuild_base else []))
    initial = work / "first/kit/compulsory-pair"
    moved = work / "neighborhood/run/kit/compulsory-pair"
    subdivisions = read(initial, "wide-finished.json")
    assert read(initial, "wide-domain.json")["grid"] == 420703
    assert subdivisions["grand_total"] == 416676
    data = read(weighted, "r7-weighted-result.json")
    lower = Fraction(data["with_seven"])
    B = Fraction(".844082046521218581802948589770767704354")
    assert data["root_profiles"] == 2534 and data["layer_pairs"] == 191069
    assert 0 < B - lower < Fraction(".000000337523")
    exclusion = work / "excluded-source"
    unpack(HERE / "excluded.zip", exclusion, "24436429d65e7b963bdfefd121ad1e280bbf839492b37e9a5f44e4dd9fc67825")
    run(exclusion, work / "excluded-proof")
    roof = read(work / "excluded-proof/previous/matched/kit/compulsory-pair", "r7-all-domain-result.json")
    relaxed = read(work / "excluded-proof/previous/kit/compulsory-pair", "r7-cap-budget.json")
    assert Fraction(roof["with_seven"]) > B and roof["root_profiles"] == 4372
    assert roof["layer_pairs"] == 401149 and roof["weighted_units"] == 1564530
    tube = read(moved, "r7-open.json")
    ctx.prec = 192
    assert arb(tube["uniform_Hausdorff_distance_lower"]) > arb(".0000000024275")
    assert tube["true_polygon"]["sizes_after_clipping"][-1] == 403
    assert arb(tube["true_polygon"]["area"]) < arb(".844082046521218581802948589770767704354")
    assert arb(relaxed["excess_at_delta"]) > arb("4.964397e-11") and relaxed["delta"] == "0.000000104"
    seed_lower = read(work / "obstruction/previous/kit/compulsory-pair", "r7-seed-free.json")["minimum_lower22"]
    result = {"fixed_cap_family": roof, "relaxed_cap_family": relaxed,
              "earlier_weighted_bound": data, "earlier_budget_gap": str(B - lower),
              "initial_host_area": read(initial, "wide-area.json")["true_area"],
              "pentagon_grid_points": 420703, "pentagon_closed_cells": subdivisions["grand_total"],
              "initial_seven_arc_residuals": [read(initial, f)["weakest_excess"] for f in
                                                ("r7-counter.json", "r7-new-counter.json", "r7-next-counter.json")],
              "repaired_five_body_polygon": read(moved, "r7-reshape-true-area.json")["r7-reshape-five.json"],
              "open_neighborhood_polygon": tube["true_polygon"], "neighborhood_radius": ".0000000024275",
              "hexagon_retaining_seed_strict_lower": seed_lower}
    assert Fraction(result["hexagon_retaining_seed_strict_lower"]) > B
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
