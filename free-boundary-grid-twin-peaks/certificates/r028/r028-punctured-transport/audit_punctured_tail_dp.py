#!/usr/bin/env python3
"""Cross-check the independent prefix DP against archived exact A2_v,B_v."""
from __future__ import annotations

import hashlib
import json
import argparse
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--r027", type=Path, required=True)
    args = ap.parse_args()
    r027 = args.r027
    out = {"status": "passed", "checks": []}
    for n in (3, 4):
        fresh_path = HERE / f"tail-n{n}.json"
        fresh = json.loads(fresh_path.read_text())
        old = json.loads((r027 / f"punctured-table-n{n}.json").read_text())
        worst_ratio = Fraction(0)
        worst_where = None
        for d, ref in zip(fresh["vertices"], old["per_vertex"]):
            assert d["A2_v_from_fibres"] == ref["A2_v"], (n, d["v"], "A")
            assert d["B_v_from_D"] == ref["B_v"], (n, d["v"], "B")
            h = d["histogram_L"]
            tails = [sum(h[i:]) for i in range(len(h))]
            for i in range(len(tails) - 1):
                if tails[i]:
                    ratio = Fraction(tails[i + 1], tails[i])
                    if ratio > worst_ratio:
                        worst_ratio = ratio
                        worst_where = {"xy": d["xy"], "t": i + 1,
                                       "T_t": tails[i], "T_t_plus_1": tails[i + 1]}
        # For odd n take the central point; for even n take the upper-left
        # member of the central 2-by-2 D4 orbit.
        cx = (n - 1) // 2
        cy = (n - 1) // 2
        centre = fresh["vertices"][cy * n + cx]
        A, B = centre["A2_v_from_fibres"], centre["B_v_from_D"]
        tail = centre["histogram_L"]
        out["checks"].append({
            "n": n,
            "archived_table_sha256": hashlib.sha256((r027 / f"punctured-table-n{n}.json").read_bytes()).hexdigest(),
            "all_vertices_match": True,
            "center_xy": centre["xy"],
            "center_B_over_A": str(Fraction(B, A)),
            "center_histogram_L": tail,
            "center_mean_L": str(Fraction(A, B)),
            "candidate_B_over_A_ge_1_5_holds_at_center": B * 5 >= A,
            "maximum_consecutive_local_tail_ratio": str(worst_ratio),
            "maximum_consecutive_local_tail_ratio_decimal": float(worst_ratio),
            "where_maximum_occurs": worst_where,
        })
    (HERE / "tail-dp-audit.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
