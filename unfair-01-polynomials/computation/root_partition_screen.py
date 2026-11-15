"""Numerically screen real root partitions of finite 0-1 polynomials.

For every C with constant and leading coefficient 1 in a stated degree range,
partition its roots into real roots and conjugate pairs.  A partition can give a
real monic factor A.  It is retained only when both factor constants are
numerically 1 and both factors have nonnegative coefficients.  Retained
non-binary cases require exact certification elsewhere; this is only a screen.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import multiprocessing as mp
import platform
import sys
import time

import numpy as np


def coeffs_from_roots(roots: list[complex]) -> np.ndarray:
    """Return ascending real coefficients of the monic root product."""
    desc = np.poly(np.asarray(roots, dtype=np.complex128))
    asc = desc[::-1]
    if np.max(np.abs(asc.imag), initial=0.0) > 2e-7:
        raise ArithmeticError("selected roots did not give a real polynomial")
    return asc.real


def root_groups(roots: np.ndarray) -> list[list[complex]] | None:
    """Group numerical roots into real singletons and conjugate pairs."""
    pending = list(complex(z) for z in roots)
    groups: list[list[complex]] = []
    while pending:
        z = pending.pop()
        if abs(z.imag) < 2e-7:
            groups.append([complex(z.real, 0.0)])
            continue
        target = z.conjugate()
        idx = min(range(len(pending)), key=lambda i: abs(pending[i] - target)) if pending else -1
        if idx < 0 or abs(pending[idx] - target) > 1e-5 * max(1.0, abs(z)):
            return None
        groups.append([z, pending.pop(idx)])
    return groups


def screen_mask(args: tuple[int, int, float]) -> dict:
    degree, mask, tol = args
    grouped = partitions = constant_one = nonnegative = fractional = 0
    candidates: list[dict] = []
    worst = {"min_coefficient": 1.0, "mask": None}
    # low-to-high internal coefficient mask; end points are fixed at one.
    c = [1.0] + [float((mask >> i) & 1) for i in range(degree - 1)] + [1.0]
    roots = np.roots(np.asarray(c[::-1]))
    groups = root_groups(roots)
    if groups is None:
        return {"root_grouped": 0, "root_partitions_tested": 0,
                "constant_one_partitions": 0, "nonnegative_partitions": 0,
                "fractional_numerical_candidates": 0, "candidates": [],
                "most_negative_constant_one_case": worst}
    grouped += 1
    g = len(groups)
    # Complementation means it suffices to test subsets containing group 0.
    for bits in range(1 << max(g - 1, 0)):
        choice = [0] + [j + 1 for j in range(g - 1) if (bits >> j) & 1]
        if len(choice) == g:
            continue
        partitions += 1
        selected = [r for j in choice for r in groups[j]]
        try:
            a = coeffs_from_roots(selected)
            b = coeffs_from_roots([r for j, group in enumerate(groups) if j not in choice for r in group])
        except ArithmeticError:
            continue
        product = np.convolve(a, b)
        if len(product) != len(c) or np.max(np.abs(product - np.asarray(c))) > 1e-5:
            continue
        if abs(a[0] - 1.0) > tol or abs(b[0] - 1.0) > tol:
            continue
        constant_one += 1
        mincoef = float(min(np.min(a), np.min(b)))
        if mincoef < worst["min_coefficient"]:
            worst = {"min_coefficient": mincoef, "mask": mask}
        if mincoef < -tol:
            continue
        nonnegative += 1
        integerish = all(abs(v - round(v)) <= tol for v in np.r_[a, b])
        if not integerish:
            fractional += 1
            candidates.append({
                "mask": mask,
                "support_C": [i for i, value in enumerate(c) if value],
                "a": [float(v) for v in a],
                "b": [float(v) for v in b],
                "min_coefficient": mincoef,
                "selection_groups": choice,
            })
    return {"root_grouped": grouped, "root_partitions_tested": partitions,
            "constant_one_partitions": constant_one, "nonnegative_partitions": nonnegative,
            "fractional_numerical_candidates": fractional, "candidates": candidates,
            "most_negative_constant_one_case": worst}


def screen_degree(degree: int, tol: float, workers: int) -> dict:
    examined = 1 << (degree - 1)
    grouped = partitions = constant_one = nonnegative = fractional = 0
    candidates: list[dict] = []
    worst = {"min_coefficient": 1.0, "mask": None}
    jobs = ((degree, mask, tol) for mask in range(examined))
    with mp.Pool(processes=workers) as pool:
        for result in pool.imap_unordered(screen_mask, jobs, chunksize=64):
            grouped += result["root_grouped"]
            partitions += result["root_partitions_tested"]
            constant_one += result["constant_one_partitions"]
            nonnegative += result["nonnegative_partitions"]
            fractional += result["fractional_numerical_candidates"]
            candidates.extend(result["candidates"])
            if result["most_negative_constant_one_case"]["min_coefficient"] < worst["min_coefficient"]:
                worst = result["most_negative_constant_one_case"]
    return {
        "degree": degree, "polynomials_examined": examined,
        "root_grouped": grouped, "root_partitions_tested": partitions,
        "constant_one_partitions": constant_one, "nonnegative_partitions": nonnegative,
        "fractional_numerical_candidates": fractional, "candidates": candidates,
        "most_negative_constant_one_case": worst,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-degree", type=int, default=4)
    parser.add_argument("--max-degree", type=int, default=14)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--workers", type=int, default=max(1, (mp.cpu_count() or 1) - 2))
    args = parser.parse_args()
    if not 2 <= args.min_degree <= args.max_degree or args.workers < 1:
        raise SystemExit("invalid degree range")
    started = time.time()
    data = {
        "algorithm": "all endpoint-one 0-1 C; numerical real root partitions",
        "certificate_status": "screen-only; no exact conclusion follows from an empty candidate list",
        "python": sys.version, "numpy": np.__version__, "platform": platform.platform(),
        "parameters": vars(args), "degrees": [],
    }
    for degree in range(args.min_degree, args.max_degree + 1):
        result = screen_degree(degree, args.tol, args.workers)
        data["degrees"].append(result)
        print(json.dumps({k: v for k, v in result.items() if k != "candidates"}, ensure_ascii=False), flush=True)
    data["elapsed_seconds"] = time.time() - started
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
