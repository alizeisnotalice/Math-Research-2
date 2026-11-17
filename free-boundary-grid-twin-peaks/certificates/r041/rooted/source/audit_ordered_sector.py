#!/usr/bin/env python3
"""Exact decimal-integer checks for ordered_sector_dp JSON output."""
import argparse, json, math, sys

P = argparse.ArgumentParser()
P.add_argument("json_file")
P.add_argument("--mode", choices=("n3", "n4"), required=True)
P.add_argument("--r040-block", help="public r040 low-branch JSON stdout (n3 gate)")
P.add_argument("--r040-image", help="public r040 direct/image JSON stdout (n3 gate)")
args = P.parse_args()
with open(args.json_file) as f: d = json.load(f)
def load_first_json(path):
    # Historical r040 stdout appends a non-JSON output_sha256 line.
    with open(path) as f: return json.JSONDecoder().raw_decode(f.read())[0]
num = lambda x: int(x)
qbs = d["ordered_q_lowb"]
A2 = num(d["A2"]); D2 = num(d["D_euclidean_squared_sum"])
by_count = sum(num(x["count"]) for x in qbs)
by_d2 = sum(num(x["D_euclidean_squared_sum"]) for x in qbs)
hist = sum(num(x["count"]) for x in d["manhattan_distance_histogram"])
ehist = sum(num(x["count"]) for x in d["euclidean_squared_distance_histogram"])
assert A2 == by_count == hist == ehist
assert D2 == by_d2
assert d["workers"] == 1 and d["hard_guard"] == "3 <= n <= 4"
qden = d["n"] * d["n"] * A2
result = {"input": args.json_file, "n": d["n"], "A1": num(d["A1"]), "A2": A2,
          "D2": D2, "Q_reduced": f"{D2//math.gcd(D2,qden)}/{qden//math.gcd(D2,qden)}",
          "hist_sum_equals_A2": True, "euclidean_bucket_sum_equals_A2": True, "table_sum_equals_A2": True,
          "table_D2_sum_equals_D2": True}
if args.mode == "n3":
    targets = {"A1": 29568, "A2": 143808, "D2": 633600,
               "rho": 15744, "lowb1_count": 79360, "lowb1_D2": 322784}
    low1 = [x for x in qbs if x["low_b"] == 1]
    got = {"A1": result["A1"], "A2": A2, "D2": D2,
           "rho": num(d["rho_manhattan_gt_3_count"]),
           "lowb1_count": sum(num(x["count"]) for x in low1),
           "lowb1_D2": sum(num(x["D_euclidean_squared_sum"]) for x in low1)}
    result["n3_gate_targets"] = targets; result["n3_gate_observed"] = got
    qnorm = result["Q_reduced"]
    result["Q_reduced"] = qnorm
    block_ok = image_ok = None
    if args.r040_block:
        ref = load_first_json(args.r040_block)
        expected = {k: v["actual_count"] for k, v in ref["rows_by_first_merge_size_and_low_branch"].items()}
        sparse = {f"q={x['q']},b={x['low_b']}": num(x["count"]) for x in qbs}
        actual = {k: sparse.get(k, 0) for k in expected}
        block_ok = actual == expected
    if args.r040_image:
        ref = load_first_json(args.r040_image)
        sparse = {f"q={x['q']},b={x['low_b']}":
                  (num(x["count"]), num(x["D_euclidean_squared_sum"]), num(x["rho_manhattan_gt_3_count"]))
                  for x in qbs}
        expected = {k: (v["direct_count"], v["direct_euclidean_squared_sum"], v["direct_Manhattan_far_count"])
                    for k, v in ref["rows"].items()}
        actual = {k: sparse.get(k, (0, 0, 0)) for k in expected}
        image_ok = (actual == expected and ref["total_euclidean_squared_sum"] == D2
                    and ref["Manhattan_far_count"] == got["rho"])
    result["r040_block_qb_counts_passed"] = block_ok
    result["r040_image_qb_count_D2_far_passed"] = image_ok
    result["n3_gate_passed"] = (got == targets and qnorm == "1100/2247" and block_ok is True and image_ok is True)
    assert result["n3_gate_passed"], (got, result["Q_reduced"])
else:
    result["legacy_global_crosscheck"] = {"A2_expected": 779299174080,
      "D2_expected": 5272309842240,
      "passed": A2 == 779299174080 and D2 == 5272309842240}
    assert result["legacy_global_crosscheck"]["passed"]
print(json.dumps(result, sort_keys=True, indent=2))
