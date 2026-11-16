#!/usr/bin/env python3
"""Exact crosschecks and rational summaries for the distance/birth joint table."""
from fractions import Fraction
import json
from pathlib import Path


BASE = Path(__file__).parent
MAIN = BASE.resolve().parents[1]
EXPECTED_A2 = {3: 143808, 4: 779299174080}
EXPECTED_DISTANCE = {
    3: {2: 63552, 3: 64512, 4: 15744},
    4: {2: 173341437440, 3: 313977112960, 4: 217620337600,
        5: 66746375360, 6: 7613910720},
}


def load(path):
    return json.loads(path.read_text())


def ints(d):
    return {int(k): int(v) for k, v in d.items()}


def summarize(data):
    n, N = data["n"], data["N"]
    joint = {int(d): ints(row) for d, row in data["joint_distance_birth"].items()}
    result = {}
    for d, row in sorted(joint.items()):
        total = sum(row.values())
        tau_sum = sum(t * c for t, c in row.items())
        gdist = {N-t: c for t, c in sorted(row.items())}
        tails = {g: sum(c for h, c in gdist.items() if h >= g)
                 for g in sorted(gdist)}
        result[str(d)] = {
            "count": str(total),
            "E_tau": str(Fraction(tau_sum, total)),
            "Pr_tau_eq_2": str(Fraction(row.get(2, 0), total)),
            "G_equals_N_minus_tau_counts": {str(g): str(c) for g, c in sorted(gdist.items())},
            "G_tail_counts_ge": {str(g): str(c) for g, c in tails.items()},
        }
    overall = ints(data["birth_counts"])
    total = sum(overall.values())
    gdist = {N-t: c for t, c in overall.items()}
    return {
        "by_distance": result,
        "overall": {
            "count": str(total),
            "E_tau": str(Fraction(sum(t*c for t, c in overall.items()), total)),
            "Pr_tau_eq_2": str(Fraction(overall.get(2, 0), total)),
            "G_equals_N_minus_tau_counts": {str(g): str(c) for g, c in sorted(gdist.items())},
            "G_tail_counts_ge": {str(g): str(sum(c for h, c in gdist.items() if h >= g))
                                    for g in sorted(gdist)},
        },
    }


def main():
    n3direct, n3dp, n4dp = (load(BASE / name) for name in
                             ("n3-direct.json", "n3-dp.json", "n4-dp.json"))
    for key in ("A2", "joint_distance_birth", "distance_counts", "birth_counts"):
        assert n3direct[key] == n3dp[key], key
    runs = []
    for data in (n3dp, n4dp):
        n = data["n"]
        assert int(data["A2"]) == EXPECTED_A2[n]
        assert ints(data["distance_counts"]) == EXPECTED_DISTANCE[n]
        assert sum(int(v) for v in data["distance_counts"].values()) == EXPECTED_A2[n]
        assert sum(int(v) for v in data["birth_counts"].values()) == EXPECTED_A2[n]
        if n == 4:
            prior = load(MAIN / "certificates/r018-root-moments/n4-d4-root-moments.json")
            assert ints(prior["distance_counts"]) == EXPECTED_DISTANCE[n]
        runs.append({"n": n, "summary": summarize(data)})

    prior_birth = load(MAIN / "certificates/r014-boundary/two-peak-birth-audit-n3.json")
    assert ints(prior_birth["second_peak_birth_step_histogram"]) == ints(n3dp["birth_counts"])
    out = {
        "status": "passed",
        "deterministic": True,
        "crosschecks": {
            "n3_all_9_factorial_direct_equals_subset_dp": True,
            "n3_existing_exact_birth_histogram_equals_joint_marginal": True,
            "n3_n4_A2_and_distance_histogram_equal_existing_certificates": True,
        },
        "runs": runs,
    }
    (BASE / "distance-birth-verified.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
