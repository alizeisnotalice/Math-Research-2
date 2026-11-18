"""Independent Python certificates for four rank-compression boundary examples."""
from collections import Counter
from fractions import Fraction
from math import gcd
import argparse
import json
from pathlib import Path


def character(v):
    return 1 - 2 * (v.bit_count() & 1)


def walsh(table):
    return {a: value for a in range(len(table))
            if (value := sum(f * character(a & x) for x, f in enumerate(table)))}


def span(points):
    vectors = {0}
    for point in points:
        vectors |= {x ^ point for x in vectors}
    return frozenset(vectors)


def all_subspaces(space):
    found = {frozenset({0})}
    new = found
    while new:
        new = {frozenset(W | {x ^ v for x in W}) for W in new for v in space} - found
        found |= new
    return sorted(found, key=lambda W: (len(W), sorted(W)))


def rank(points):
    return len(span(points)).bit_length() - 1


def verify(table, core):
    n = len(table).bit_length() - 1
    ambient = range(1 << n)
    original = walsh(table)
    divisor = gcd(*original.values())
    odd = {a for a, value in original.items() if value // divisor & 1}
    offset = min(core)
    direction = frozenset(a ^ offset for a in core)
    assert span(direction) == direction and set(core) <= original.keys()
    signs = [f * character(x & offset) for x, f in enumerate(table)]
    spectrum = walsh(signs)
    quotient = {point: min(point ^ w for w in direction) for point in ambient}
    bands = {gamma: sorted(a for a in spectrum if quotient[a] == gamma)
             for gamma in set(quotient.values()) if gamma}
    bands = {gamma: row for gamma, row in bands.items() if row}
    external_rank = rank(bands)
    q = sum(map(len, bands.values()))
    basis = []
    for gamma in sorted(bands):
        if gamma not in span(basis):
            basis.append(gamma)
    differences = {a ^ bands[gamma][0] for gamma in basis for a in bands[gamma][1:]}
    budget = sum(len(bands[gamma]) - 1 for gamma in basis)
    M = len(direction)
    profiles = {}
    for W in all_subspaces(direction):
        L = M // len(W)
        normal = frozenset(x for x in ambient if all(character(x & w) == 1 for w in W))
        shifts = sorted({min(x ^ u for u in normal) for x in ambient})
        frequency = sorted({min(a ^ w for w in W) for a in ambient})
        p = profiles.setdefault(L, {"ranks": Counter(), "least_sparsity_for_full_rank": None,
                                    "collision_free_kernels": 0, "full_rank_witness": None})
        p["collision_free_kernels"] += not W & differences
        for x in shifts:
            restricted = {}
            for a in frequency:
                value = sum(signs[x ^ u] * character(a & u) for u in normal)
                assert value * len(W) == sum(spectrum.get(a ^ w, 0) * character((a ^ w) & x) for w in W)
                if value:
                    restricted[a] = value
            kept = rank(quotient[a] for a in restricted)
            assert kept <= external_rank and len(restricted) <= q + L
            if not W & differences:
                assert kept == external_rank
            p["ranks"][kept] += 1
            if kept == external_rank and (p["least_sparsity_for_full_rank"] is None
                                           or len(restricted) < p["least_sparsity_for_full_rank"]):
                p["least_sparsity_for_full_rank"] = len(restricted)
                p["full_rank_witness"] = {"kernel": sorted(W), "physical_shift": x,
                                          "restricted_walsh": restricted}
    for L, profile in profiles.items():
        histogram = profile["ranks"]
        loss = Fraction(sum((external_rank - kept) * count for kept, count in histogram.items()),
                        sum(histogram.values()))
        if M > 1:
            cost = Fraction(budget * (M // L - 1), M - 1)
            assert 2 * loss <= cost
            assert external_rank - loss >= Fraction(external_rank ** 2, external_rank + cost)
        profile["mean_rank_loss"] = str(loss)
    return {"n": n, "negative_points": [i for i, f in enumerate(table) if f < 0],
            "original_walsh": original, "actual_highest_odd_layer": sorted(odd),
            "core": sorted(core), "M": M, "q": q, "external_rank": external_rank,
            "basis": basis, "selected_bands": {gamma: bands[gamma] for gamma in basis},
            "difference_set": sorted(differences), "multiplicity_budget": budget,
            "profiles_by_L": profiles}


def run():
    small = verify([1 - 2 * ((23 >> x) & 1) for x in range(16)], {1, 9})
    assert small["profiles_by_L"][1]["ranks"] == {1: 1, 2: 1}
    dedup = verify([1 - 2 * (x == 0) for x in range(16)], {0, 1})
    assert dedup["profiles_by_L"][1]["mean_rank_loss"] == "3/2"
    assert len(dedup["difference_set"]) == 1

    # Address variables occupy the low bits, data use the high bits.
    two_selector = []
    for x in range(64):
        s, y, z = x & 3, (x >> 2) & 3, x >> 4
        local = (1 - 2 * (y == 0), character(z & 1), character(z & 2), 1)
        two_selector.append(local[s])
    highest_loss = verify(two_selector, set(range(16)))
    assert highest_loss["actual_highest_odd_layer"] == list(range(16))
    assert max(highest_loss["profiles_by_L"][1]["ranks"]) == 1
    assert max(highest_loss["profiles_by_L"][2]["ranks"]) == 2

    one_selector = []
    for x in range(128):
        t, y, z, v = x & 1, (x >> 1) & 7, (x >> 4) & 3, x >> 6
        one_selector.append(character(v ^ (z == 3)) if t else 1 - 2 * (y == 0))
    highest_dedup = verify(one_selector, set(range(16)))
    assert highest_dedup["actual_highest_odd_layer"] == list(range(16))
    assert highest_dedup["profiles_by_L"][1]["mean_rank_loss"] == "3/2"
    assert len(highest_dedup["difference_set"]) == 1

    return {"smallest_truth_mask_obstruction": small, "smallest_deduplicated_bound_failure": dedup,
            "highest_layer_full_collapse_rank_loss": highest_loss,
            "highest_layer_deduplicated_bound_failure": highest_dedup}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    results = run()
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    for name, result in results.items():
        print(name, {L: row["mean_rank_loss"] for L, row in result["profiles_by_L"].items()})
