#!/usr/bin/env python3
"""Independent integer checks for beta_zero_layered archived JSON outputs.

This script deliberately does not execute a large enumeration.  It verifies
the stored n=7 result's orbit coverage and arithmetic aggregation, compares
n=3/5/6 totals with independently established values, and emits a small JSON
certificate.  It uses only the Python standard library and exact integers.
"""
import argparse
import json
from pathlib import Path


EXPECTED_B0 = {
    3: 64448,
    5: 5709646459857341216640,
    6: 275686968906490917891735764030554880,
    7: 131529421404957780122938007473728844832061946591752192,
}
A2_7 = 137077928808481869759712347602605837601908319657172480


def load(path: Path):
    with path.open() as fh:
        return json.load(fh)


def check_small(base: Path):
    checked = {}
    for n in (3, 5, 6):
        x = load(base / f"beta-zero-layered-n{n}.json")
        assert x["n"] == n
        assert int(x["beta0_total"]) == EXPECTED_B0[n]
        assert x["directed_edge_coverage"] == 4 * n * (n - 1)
        checked[str(n)] = x["beta0_total"]
    return checked


def check_n7(base: Path):
    x = load(base / "beta-zero-layered-n7.json")
    assert x["n"] == 7 and len(x["edge_orbits"]) == 12
    assert x["undirected_edge_coverage"] == 84
    assert x["directed_edge_coverage"] == 168 == 4 * 7 * 6
    assert sum(2 * o["multiplicity"] for o in x["edge_orbits"]) == 168
    total = sum(2 * o["multiplicity"] * int(o["beta0_one_direction"])
                for o in x["edge_orbits"])
    assert total == int(x["beta0_total"]) == EXPECTED_B0[7]
    support = A2_7 - total
    assert support > 0
    return {"beta0": str(total), "A2": str(A2_7), "beta_positive_support": str(support)}


def check_n7_orbit_replay(orbit_dir: Path):
    manifest = load(orbit_dir / "manifest.json")
    assert manifest["n"] == 7 and manifest["orbits_requested"] == 12
    assert manifest["status"] == "all_orbits_completed"
    assert manifest["source_sha256_before"] == manifest["source_sha256_after"]
    assert manifest["binary_sha256_before"] == manifest["binary_sha256_after"]
    assert len(manifest["runs"]) == 12
    total = 0
    coverage = 0
    for expected, run in enumerate(manifest["runs"]):
        assert run["orbit"] == expected and run["exit_code"] == 0
        assert run["post_exit_json_check"] is True
        x = load(orbit_dir / run["json"])
        assert x["selected_orbit"] == expected and len(x["edge_orbits"]) == 1
        item = x["edge_orbits"][0]
        assert x["directed_edge_coverage"] == 2 * item["multiplicity"]
        total += 2 * item["multiplicity"] * int(item["beta0_one_direction"])
        coverage += x["directed_edge_coverage"]
        timing = (orbit_dir / run["time"]).read_text()
        assert "maximum resident set size" in timing and "peak memory footprint" in timing
    assert coverage == 168 and total == EXPECTED_B0[7]
    return {"orbits": 12, "directed_coverage": coverage,
            "beta0": str(total),
            "source_sha256": manifest["source_sha256_before"],
            "binary_sha256": manifest["binary_sha256_before"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, default=Path(__file__).parent)
    ap.add_argument("--certificate", type=Path, required=True)
    ap.add_argument("--orbit-dir", type=Path)
    args = ap.parse_args()
    out = {"small_exact_crosschecks": check_small(args.base),
           "n7_exact_checks": check_n7(args.base),
           "status": "passed"}
    if args.orbit_dir:
        out["n7_orbit_replay"] = check_n7_orbit_replay(args.orbit_dir)
    args.certificate.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
