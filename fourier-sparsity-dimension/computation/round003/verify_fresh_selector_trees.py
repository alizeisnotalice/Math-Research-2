#!/usr/bin/env python3
"""Exact finite audit for A28 fresh-selector binary trees."""
import argparse
import hashlib
import itertools
import json
import platform
import random
import sys
import time

def popcount(value):
    return bin(value).count("1")

def rank(vectors):
    pivots = {}
    for vector in vectors:
        x = vector
        while x:
            bit = x.bit_length() - 1
            if bit in pivots:
                x ^= pivots[bit]
            else:
                pivots[bit] = x
                break
    return len(pivots)

def fwht(values):
    out = values[:]
    step = 1
    while step < len(out):
        for base in range(0, len(out), 2 * step):
            for j in range(base, base + step):
                left, right = out[j], out[j + step]
                out[j], out[j + step] = left + right, left - right
        step *= 2
    return out

def spectrum(values):
    return {frequency for frequency, coefficient in enumerate(fwht(values)) if coefficient}

def lift_rank(frequencies, dimension):
    return rank([frequency | (1 << dimension) for frequency in frequencies])

def make_leaf(label, nselectors, p, ambient_dimension):
    b = label & ((1 << p) - 1)
    sign = 1 if label < (1 << p) else -1
    values = []
    for x in range(1 << ambient_dimension):
        y = x >> nselectors
        values.append(sign if popcount(b & y) % 2 == 0 else -sign)
    freqs = spectrum(values)
    return {"values": values, "freqs": freqs, "a": lift_rank(freqs, ambient_dimension)}

def selector_free(stat, nselectors):
    selector_mask = (1 << nselectors) - 1
    return all((frequency & selector_mask) == 0 for frequency in stat["freqs"])

def make_node(selector, left, right, nselectors, ambient_dimension, path, audit):
    values = [left["values"][x] if ((x >> selector) & 1) == 0 else right["values"][x]
              for x in range(1 << ambient_dimension)]
    freqs = spectrum(values)
    actual_a = lift_rank(freqs, ambient_dimension)
    h = left["a"] + right["a"] - rank(
        [frequency | (1 << ambient_dimension) for frequency in left["freqs"]] +
        [frequency | (1 << ambient_dimension) for frequency in right["freqs"]])
    predicted_a = left["a"] + right["a"] - h + 1
    # Fresh left/right selector sets are disjoint: an intersection can only be pure data.
    common_pure_data = len(left["freqs"] & right["freqs"])
    loss = 2 * len(left["freqs"]) + 2 * len(right["freqs"]) - len(freqs)
    if not 2 * common_pure_data <= loss <= 3 * common_pure_data:
        raise AssertionError(("sparsity recurrence", path, common_pure_data, loss))
    exceptional = selector_free(left, nselectors) and selector_free(right, nselectors)
    if not exceptional and actual_a != predicted_a:
        raise AssertionError(("affine recurrence", path, actual_a, predicted_a, h))
    if not exceptional and actual_a == left["a"] + right["a"] - h:
        raise AssertionError(("altered coefficient accepted", path))
    affine_dimension = actual_a - 1
    if 4 * len(freqs) < (affine_dimension + 2) ** 2:
        raise AssertionError(("quadratic bound", path, len(freqs), affine_dimension))
    audit["nodes"] += 1
    audit["loss_histogram"][str(loss)] = audit["loss_histogram"].get(str(loss), 0) + 1
    audit["exception_counts"]["selectorfree_pair" if exceptional else "nonexception"] += 1
    audit["minimum_slack"] = min(audit["minimum_slack"], 4 * len(freqs) - (affine_dimension + 2) ** 2)
    audit["altered_coefficient_rejected"] = audit["altered_coefficient_rejected"] or not exceptional
    return {"values": values, "freqs": freqs, "a": actual_a}

def fresh_audit():
    return {"nodes": 0, "loss_histogram": {}, "exception_counts": {"selectorfree_pair": 0, "nonexception": 0},
            "minimum_slack": 10 ** 9, "altered_coefficient_rejected": False}

def build(shape, labels, nselectors, p, audit):
    ambient_dimension = nselectors + p
    iterator = iter(labels)
    def walk(tree, path):
        if tree == "leaf":
            return make_leaf(next(iterator), nselectors, p, ambient_dimension)
        selector, left_tree, right_tree = tree
        left = walk(left_tree, path + "L")
        right = walk(right_tree, path + "R")
        return make_node(selector, left, right, nselectors, ambient_dimension, path, audit)
    root = walk(shape, "")
    if next(iterator, None) is not None:
        raise ValueError("extra labels")
    return root

def run_family(name, shape, labels_iter, nselectors, p):
    audit = fresh_audit()
    roots = 0
    best = None
    for labels in labels_iter:
        root = build(shape, labels, nselectors, p, audit)
        candidate = (len(root["freqs"]), root["a"] - 1, list(labels))
        if best is None or candidate[0] * (best[1] + 2) ** 2 < best[0] * (candidate[1] + 2) ** 2:
            best = candidate
        roots += 1
    return {"name": name, "roots": roots, "root_best": {"s": best[0], "affdim": best[1], "labels": best[2]}, "audit": audit}

def edge_cases():
    shape = (0, "leaf", "leaf")
    # label = b + 4*negative-sign, p=2
    cases = {"equal_constant": (0, 0), "opposite_constant_root_shift": (0, 4),
             "equal_character": (3, 3), "opposite_character_root_shift": (3, 7)}
    answer = {}
    for name, labels in cases.items():
        audit = fresh_audit()
        root = build(shape, labels, 1, 2, audit)
        answer[name] = {"s": len(root["freqs"]), "affdim": root["a"] - 1, "audit": audit}
    return answer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--depth3-samples", type=int, default=2000)
    args = parser.parse_args()
    begun = time.monotonic()
    balanced4 = (0, (1, "leaf", "leaf"), (2, "leaf", "leaf"))
    comb4 = (0, "leaf", (1, "leaf", (2, "leaf", "leaf")))
    result = {"method": "integer_truth_table_fwht_and_f2_elimination", "seed": 20261116,
              "four_leaf_balanced": run_family("balanced4", balanced4, itertools.product(range(8), repeat=4), 3, 2),
              "four_leaf_comb": run_family("comb4", comb4, itertools.product(range(8), repeat=4), 3, 2)}
    balanced8 = (0, (1, (2, "leaf", "leaf"), (3, "leaf", "leaf")),
                 (4, (5, "leaf", "leaf"), (6, "leaf", "leaf")))
    generator = random.Random(20261116)
    labels = [tuple(generator.randrange(16) for _ in range(8)) for _ in range(args.depth3_samples)]
    result["eight_leaf_balanced_sample"] = run_family("balanced8", balanced8, labels, 7, 3)
    result["edge_cases"] = edge_cases()
    result["elapsed_seconds"] = time.monotonic() - begun
    result["python"] = sys.version
    result["platform"] = platform.platform()
    for group in (result["four_leaf_balanced"], result["four_leaf_comb"], result["eight_leaf_balanced_sample"]):
        if not group["audit"]["altered_coefficient_rejected"]:
            raise AssertionError("altered coefficient test was vacuous")
    blob = json.dumps(result, indent=2, sort_keys=True) + "\n"
    open(args.out, "w").write(blob)
    print(json.dumps({"status": "passed", "sha256": hashlib.sha256(blob.encode()).hexdigest(), "elapsed_seconds": result["elapsed_seconds"]}))

if __name__ == "__main__":
    main()
