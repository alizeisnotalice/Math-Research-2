#!/usr/bin/env python3
"""Exact affine half-window codes for archived, tree, and four-variable flat spectra."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import hashlib
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import platform
import time


def parity(n):
    return n.bit_count() & 1


def indexes(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def fwht(input):
    data = list(input)
    w = 1
    while w < len(data):
        for start in range(0, len(data), 2 * w):
            for i in range(start, start + w):
                a, b = data[i], data[i + w]
                data[i], data[i + w] = a + b, a - b
        w *= 2
    return data


def independent(words):
    pivots = {}
    basis = []
    for x in words:
        y = x
        while y:
            k = y.bit_length() - 1
            if k in pivots:
                y ^= pivots[k]
            else:
                pivots[k] = y
                basis.append(x)
                break
    return basis


def solve(words, target):
    pivots = {}
    for x, t in zip(words, target):
        while x:
            p = x.bit_length() - 1
            if p not in pivots:
                pivots[p] = x, t
                break
            row, bit = pivots[p]
            x ^= row
            t ^= bit
        if not x:
            assert not t
    y = 0
    for p, (row, bit) in sorted(pivots.items()):
        y |= (bit ^ parity(row & y)) << p
    assert all(parity(x & y) == t for x, t in zip(words, target))
    return y


def span(rows):
    out = [0]
    for row in rows:
        out += [w ^ row for w in out]
    return out


def packed_evaluations(points, n):
    aug = [p | (1 << n) for p in points]
    rows = [sum(((p >> bit) & 1) << i for i, p in enumerate(aug)) for bit in range(n + 1)]
    bases = independent(rows)
    reps = [1 << rows.index(b) for b in bases]
    return span(bases), span(reps)


@cache
def physical_masks(n):
    full = (1 << (1 << n)) - 1
    sides = [full // ((1 << (1 << j)) + 1) for j in range(n)]
    return full, sides


def character_negatives(frequency, n):
    full, sides = physical_masks(n)
    out = 0
    for j in indexes(frequency):
        out ^= full ^ sides[j]
    return out


def translate_truth(mask, q, n):
    _, sides = physical_masks(n)
    for j in indexes(q):
        shift = 1 << j
        lower = sides[j]
        mask = ((mask & lower) << shift) | ((mask >> shift) & lower)
    return mask


def minimum_checks(points, signs, n, minimum, word_reps, truth=None, physical=None):
    dim = minimum.bit_length() - 1
    assert minimum == 1 << dim
    full, _ = physical_masks(n) if truth is not None else (None, None)
    for mask, rep in word_reps:
        locations = list(indexes(mask))
        target = [points[i] for i in locations]
        origin = target[0]
        basis = independent(x ^ origin for x in target)
        assert len(basis) == dim and {x ^ origin for x in span(basis)} == set(target)
        sign_t = [(1 - signs[i]) // 2 for i in locations]
        z = solve([x | (1 << n) for x in target], sign_t)
        x0, c = z & ((1 << n) - 1), z >> n
        q, b = rep & ((1 << n) - 1), rep >> n
        assert all(parity(q & x) ^ b for x in target)
        if truth is not None:
            actual = truth ^ translate_truth(truth, q, n) ^ (full if b else 0)
            support = full
            for a in basis:
                half = character_negatives(a, n)
                support &= half if parity(a & x0) else full ^ half
            assert actual == support
            assert actual.bit_count() * minimum == 1 << n
            sign = character_negatives(origin, n) ^ (full if c ^ parity(origin & x0) else 0)
            assert not (truth ^ sign) & actual
        elif physical is not None:
            for point in {x0, x0 ^ q} | {x0 ^ (1 << i) for i in range(n)}:
                value = sum(s * (1 - 2 * parity(a & point)) for a, s in zip(target, [signs[i] for i in locations]))
                assert 2 * value == minimum * (physical(point) - (1 - 2 * b) * physical(point ^ q))


@cache
def rref(n, r):
    out = []
    for p in combinations(range(n), r):
        free = [(i, j) for i, bit in enumerate(p) for j in range(bit + 1, n) if j not in p]
        for mask in range(1 << len(free)):
            rows = [1 << j for j in p]
            for k, (i, j) in enumerate(free):
                rows[i] |= ((mask >> k) & 1) << j
            out.append(rows)
    return out


@cache
def small_exact(words, d):
    best, witness = max(w.bit_length() for w in words) + 1, []
    for generators in rref(len(words).bit_length() - 1, d):
        w = [words[q] for q in generators]
        u = 0
        for x in w:
            u |= x
        if u.bit_count() < best:
            best, witness = u.bit_count(), w
    return {"lower": best, "upper": best, "basis_masks": [hex(w) for w in witness], "proof": "all RREF subcodes"}


def shell_bounds(minimum_words, w1, w2, dimensions, extra=None):
    found = {}
    minset = set(minimum_words)
    prev_max = 1
    for d in (2, 3):
        if dimensions < d:
            continue
        best = 1 << 30
        max_minimum = prev_max
        witness = ()
        for cand in combinations(minimum_words, d):
            codes = span(cand)
            if 0 in codes[1:]:
                continue
            max_minimum = max(max_minimum, len(minset.intersection(codes)))
            merged = 0
            for c in cand:
                merged |= c
            if merged.bit_count() < best:
                best, witness = merged.bit_count(), cand
        if not witness and extra is not None:
            witness = tuple(independent(extra)[:d])
            total = 0
            for x in witness:
                total |= x
            best = total.bit_count()
        assert witness and len(witness) == d
        slots = 1 << (d - 1)
        lower = (max_minimum * w1 + ((1 << d) - 1 - max_minimum) * w2 + slots - 1) // slots
        assert lower <= best
        found[str(d)] = {"lower": lower, "upper": best, "maximum_minimum_words_in_subcode": max_minimum,
                         "basis_masks": [hex(x) for x in witness], "proof": "minimum-shell counting and weight averaging"}
        prev_max = max_minimum
    return found


def describe(code_hist, word_reps, points, signs, n, dimension, norm, *, truth=None, actual=None, all_words=None):
    lengths = sorted(k for k in code_hist if k)
    w1, w2 = (lengths + [None])[:2]
    assert w1 >= norm and len(points) == norm * norm
    assert sum(code_hist.values()) == 1 << dimension
    assert all(code_hist[k] == code_hist[len(points) - k] for k in code_hist)
    assert len(word_reps) == code_hist[w1] and len({w for w, _ in word_reps}) == len(word_reps)
    if w1 == norm:
        minimum_checks(points, signs, n, w1, word_reps, truth, actual)
    mins = [w for w, _ in word_reps]
    counts, min_dimensions = Counter(), Counter()
    for u in mins:
        vectors = [points[i] for i in indexes(u)]
        dim = len(independent(x ^ vectors[0] for x in vectors))
        min_dimensions[str(dim)] += 1
    for u, v in combinations(mins, 2):
        both = u & v
        if not both:
            counts["empty"] += 1
            continue
        target = [points[i] for i in indexes(both)]
        dim = len(independent(x ^ target[0] for x in target))
        if w1 == norm:
            assert len(target) == (1 << dim)
        counts[str(dim)] += 1
    if dimension <= 7:
        assert all_words is not None
        bounds = {str(d): small_exact(tuple(all_words), d) for d in (2, 3) if dimension >= d}
    else:
        bounds = shell_bounds(mins, w1, w2, dimension, all_words)
    return {
        "flat_N": norm, "s": len(points), "code_dimension": dimension,
        "nonzero_minimum_weight": w1, "next_positive_weight": w2,
        "minimum_words": code_hist[w1], "next_words": code_hist[w2] if w2 is not None else 0,
        "weight_histogram": dict(sorted(code_hist.items())), "minimum_support_affine_dims": dict(min_dimensions),
        "minimum_pair_intersection_affine_dims": dict(counts),
        "number_of_N_equality_signchecks": len(word_reps) if w1 == norm else 0,
        "physical_derivative_check": "all points by bit permutation" if truth is not None else "explicit tree values at witness probes",
        "generalized_weights": bounds,
    }


def full_code(points, signs, n, norm, truth):
    words, reps = packed_evaluations(points, n)
    hist = Counter(w.bit_count() for w in words)
    w1 = min(w for w in hist if w)
    minima = [(w, q) for w, q in zip(words, reps) if w.bit_count() == w1]
    return describe(hist, minima, points, signs, n, len(words).bit_length() - 1, norm,
                    truth=truth, all_words=words)


def graph11(row):
    points, signs = row["support"], row["coefficient_numerators_over_8"]
    spectrum = [0] * (1 << 15)
    for a, c in zip(points, signs):
        spectrum[a] = c
    truth = fwht(spectrum)
    assert {abs(y) for y in truth} == {8}
    return {"period_dimension": row["period_dimension"], "archived_row_index": row["row_index"],
            **full_code(points, signs, 15, 8, sum((v < 0) << i for i, v in enumerate(truth)))}


LINES = ((0, 1), (2, 4), (3, 7), (5, 6))


def tree49(height):
    B, block = 1 << (height - 1), (1 << height) - 1
    N, n = 4 * B, 3 + 4 * block
    locals = []
    global_spectrum = []
    for i, (a, b) in enumerate(LINES):
        H = [q for q in range(8) if not parity(q & (a ^ b))]
        group = []
        for terminal in range(B - 1, block):
            variants = [(1 << terminal, 0)]
            v = terminal
            while v:
                parent = (v - 1) // 2
                bit = int(v == 2 * parent + 2)
                variants += [(q ^ (1 << parent), s ^ bit) for q, s in variants]
                v = parent
            for point, s in variants:
                for addr in H:
                    group.append((point, addr))
                    global_spectrum.append((addr | (point << (3 + block * i)), 1 - 2 * (s ^ parity(addr & a))))
        locals.append(group)
    points, signs = zip(*global_spectrum)
    assert len(points) == len(set(points)) == N * N and len(independent(points)) == n

    def evaluate(x):
        i = next(i for i, line in enumerate(LINES) if x % 8 in line)
        v, local = 0, x >> (3 + block * i)
        while v < B - 1:
            v = 2 * v + 1 + ((local >> v) & 1)
        return 1 - 2 * ((local >> v) & 1)

    if height <= 2:
        exact = fwht(evaluate(x) for x in range(1 << n))
        reference = dict(global_spectrum)
        assert all(val * N == (1 << n) * reference.get(i, 0) for i, val in enumerate(exact))
        return full_code(points, signs, n, N, sum((evaluate(x) < 0) << x for x in range(1 << n)))

    plain = [
        span([sum(((local >> j) & 1) << k for k, (local, _) in enumerate(g)) for j in range(block)])
        for g in locals
    ]
    assert all(len(set(words)) == 1 << block for words in plain)
    histogram = Counter()
    minima = []
    group_size = len(locals[0])
    for q in range(8):
        masks = [sum(parity(q & a) << j for j, (_, a) in enumerate(g)) for g in locals]
        grouped = []
        combined = Counter({0: 1})
        for i in range(4):
            smalls = []
            weights = Counter()
            for local_q, word in enumerate(plain[i]):
                word ^= masks[i]
                k = word.bit_count()
                weights[k] += 1
                if k <= N:
                    smalls.append((k, local_q, word))
            grouped.append(smalls)
            next_counts = Counter()
            for a, count in combined.items():
                for b, ways in weights.items():
                    next_counts[a + b] += count * ways
            combined = next_counts
        histogram.update(combined)

        def collect(i, rem, representative, word):
            if i == 4:
                if rem == 0:
                    minima.append((word, representative))
                return
            for k, local_q, local_word in grouped[i]:
                if k <= rem:
                    collect(i + 1, rem - k, representative ^ (local_q << (3 + block * i)),
                            word | (local_word << (i * group_size)))
        collect(0, N, q, 0)
    assert min(k for k in histogram if k) == N and histogram[N] == len(minima)
    for x in [0, (1 << n) - 1] + [1 << j for j in range(n)]:
        assert sum(c * (1 - 2 * parity(a & x)) for a, c in global_spectrum) == N * evaluate(x)
    return describe(histogram, minima, points, signs, n, n, N, actual=evaluate)


def four_shard(bounds):
    lo, hi = bounds
    classes, witnesses = Counter(), {}
    checked = Counter()
    for mask in range(lo, hi):
        W = fwht(1 - 2 * ((mask >> i) & 1) for i in range(16))
        magnitudes = {abs(v) for v in W if v}
        if len(magnitudes) != 1:
            continue
        points = [i for i, v in enumerate(W) if v]
        N = isqrt(len(points))
        assert 16 == N * next(iter(magnitudes))
        signs = [1 - 2 * (W[x] < 0) for x in points]
        report = full_code(points, signs, 4, N, mask)
        key = (report["s"], report["code_dimension"], report["nonzero_minimum_weight"], report["next_positive_weight"])
        classes[key] += 1
        witnesses.setdefault(key, {"truth_mask": hex(mask), **report})
        checked["distinct_affine_codewords_across_flat_functions"] += 1 << report["code_dimension"]
        checked["equality_windows_physically_verified"] += report["number_of_N_equality_signchecks"]
    return classes, witnesses, checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round011", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    start = time.perf_counter()
    graph_rows = json.loads(args.round011.read_bytes())["rows"]
    assert len(graph_rows) == 20
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        graph = list(pool.map(graph11, graph_rows))
        tree = list(pool.map(tree49, range(1, 5)))
        classes, witnesses, stats = Counter(), {}, Counter()
        for part, example, counts in pool.map(four_shard, [(i, i + 4096) for i in range(0, 65536, 4096)]):
            classes.update(part)
            stats.update(counts)
            for key, record in example.items():
                witnesses.setdefault(key, record)
    result = {
        "python": platform.python_version(), "workers": args.workers,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "round011_input_sha256": hashlib.sha256(args.round011.read_bytes()).hexdigest(),
        "graph_64_point_spectra": graph,
        "r49_four_line_spectra_by_height": {i: x for i, x in enumerate(tree, 1)},
        "r49_counting_method": "H<=2 direct evaluation code, H=3,4 full compressed enumeration by fixed address and four independent block codes",
        "four_variable_truth_tables": 65536, "four_variable_flat_functions": sum(classes.values()),
        "four_variable_classes": [dict(function_count=classes[k], **witnesses[k]) for k in sorted(classes)],
        "four_variable_checks": dict(stats), "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
