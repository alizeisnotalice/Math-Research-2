#!/usr/bin/env python3
"""CPU-only integer matrix cache builder for positive_certificate_search.

The sparse builder calls ``Poly(...).as_dict()`` once per column and scatters
the resulting exponent dictionary into the original monomial-basis order.  It
does not call an LP solver and never creates a torch/MPS object.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
from itertools import combinations, product
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as S


HERE = Path(__file__).resolve().parent
FROZEN = HERE / "positive_certificate_search_frozen.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def load_helpers():
    spec = importlib.util.spec_from_file_location("pcs_frozen", FROZEN)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_support(text: str) -> list[int]:
    return [int(value) for value in text.split(",")]


def structure(pcs, I: list[int], J: list[int], multiplier_degree: int, atom_degree: int):
    """Reproduce source column and basis order exactly, before matrix filling."""
    pcs.validate_supports(I, J)
    D = pcs.diagonals(I, J)
    units = pcs.derived_units(I, J, D)
    keys = [("a", x) for x in I] + [("b", x) for x in J]
    strict_keys = [key for key in keys if key not in units]
    symbol_result = S.symbols(" ".join(pcs.variable_name(*key) for key in strict_keys)) if strict_keys else ()
    symbols = ((symbol_result,) if len(strict_keys) == 1 else tuple(symbol_result))
    values = {key: S.Integer(1) for key in units}
    values.update({key: symbol for key, symbol in zip(strict_keys, symbols)})
    equations = {k: S.expand(sum(values[("a", i)] * values[("b", j)] for i, j in pairs) - 1)
                 for k, pairs in D.items()}
    equations = {k: value for k, value in equations.items() if value != 0}
    multipliers = pcs.monomials(symbols, multiplier_degree)
    product_columns = [(k, monomial, S.expand(monomial * equation))
                       for k, equation in equations.items() for monomial in multipliers]
    atoms = []
    for size in range(1, atom_degree + 1):
        for picks in combinations(range(len(symbols)), size):
            for signs in product((0, 1), repeat=size):
                positive = tuple(str(symbols[index]) for index, sign in zip(picks, signs) if sign == 0)
                one_minus = tuple(str(symbols[index]) for index, sign in zip(picks, signs) if sign == 1)
                atom = S.expand(S.prod(symbols[index] if sign == 0 else 1 - symbols[index]
                                       for index, sign in zip(picks, signs)))
                atoms.append((positive, one_minus, atom))
    max_degree = max(S.Poly(poly, *symbols).total_degree()
                     for _, _, poly in product_columns + [(None, None, atom) for _, _, atom in atoms])
    basis = pcs.monomials(symbols, max_degree)
    return {"D": D, "units": units, "symbols": symbols, "equations": equations,
            "product_columns": product_columns, "atoms": atoms, "basis": basis,
            "max_degree": max_degree}


def sparse_columns(columns, symbols, basis) -> np.ndarray:
    """One Poly.as_dict per column; basis order is the source monomial order."""
    row_of = {tuple(monomial.as_powers_dict().get(symbol, S.Integer(0)) for symbol in symbols): row
              for row, monomial in enumerate(basis)}
    result = np.zeros((len(basis), len(columns)), dtype=np.int64)
    for column, (_, _, expression) in enumerate(columns):
        terms = S.Poly(expression, *symbols).as_dict()
        for exponents, coefficient in terms.items():
            assert coefficient.is_Integer, (column, exponents, coefficient)
            result[row_of[exponents], column] = int(coefficient)
    return result


def old_columns(pcs, columns, symbols, basis) -> np.ndarray:
    """Intentionally slow source-equivalent per-monomial coefficient query."""
    return np.asarray([[int(value) for value in pcs.polynomial_vector(expression, symbols, basis)]
                       for _, _, expression in columns], dtype=np.int64).T


def layout_metadata(pcs, data, I, J, multiplier_degree, atom_degree) -> dict[str, object]:
    product_layout = [{"diagonal": diagonal, "monomial": pcs.encode_monomial(monomial)}
                      for diagonal, monomial, _ in data["product_columns"]]
    atom_layout = [{"positive": list(positive), "one_minus": list(one_minus)}
                   for positive, one_minus, _ in data["atoms"]]
    basis_layout = [pcs.encode_monomial(monomial) for monomial in data["basis"]]
    return {
        "I": I, "J": J, "multiplier_degree": multiplier_degree, "atom_degree": atom_degree,
        "strict_variables": [str(symbol) for symbol in data["symbols"]],
        "unit_variables": sorted(pcs.variable_name(*key) for key in data["units"]),
        "equation_diagonals": list(data["equations"]),
        "max_total_degree": data["max_degree"],
        "basis_hash": canonical_hash(basis_layout),
        "product_column_layout_hash": canonical_hash(product_layout),
        "atom_column_layout_hash": canonical_hash(atom_layout),
    }


def build(I, J, multiplier_degree, atom_degree, compare_old: bool):
    pcs = load_helpers()
    started = perf_counter()
    data = structure(pcs, I, J, multiplier_degree, atom_degree)
    if not data["symbols"]:
        raise ValueError("no strict variables: no nontrivial matrix cache")
    L = sparse_columns(data["product_columns"], data["symbols"], data["basis"])
    R = sparse_columns(data["atoms"], data["symbols"], data["basis"])
    matrix = np.hstack((L, -L, -R))
    normalization = np.concatenate((np.zeros(2 * L.shape[1], dtype=np.int64), np.ones(R.shape[1], dtype=np.int64)))
    Aeq = np.vstack((matrix, normalization))
    elapsed = perf_counter() - started
    old_check = {"performed": False}
    if compare_old:
        old_started = perf_counter()
        old_L = old_columns(pcs, data["product_columns"], data["symbols"], data["basis"])
        old_R = old_columns(pcs, data["atoms"], data["symbols"], data["basis"])
        old_Aeq = np.vstack((np.hstack((old_L, -old_L, -old_R)), normalization))
        assert np.array_equal(Aeq, old_Aeq), "sparse matrix differs from old coefficient-query build"
        old_check = {"performed": True, "all_entries_equal": True,
                     "old_build_seconds": perf_counter() - old_started,
                     "old_integer_matrix_sha256": hashlib.sha256(old_Aeq.tobytes()).hexdigest()}
    metadata = layout_metadata(pcs, data, I, J, multiplier_degree, atom_degree)
    metadata.update({
        "matrix_shape": [int(Aeq.shape[0]), int(Aeq.shape[1])],
        "polynomial_rows": int(matrix.shape[0]), "normalization_rows": 1,
        "free_multiplier_columns": int(L.shape[1]), "positive_atom_columns": int(R.shape[1]),
        "integer_matrix_sha256": hashlib.sha256(Aeq.tobytes()).hexdigest(),
        "sparse_build_seconds": elapsed, "old_matrix_comparison": old_check,
    })
    return Aeq, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--I", default="0,1,2,4,6,9,11,16")
    parser.add_argument("--J", default="0,1,3,4,6,8,11,13,18")
    parser.add_argument("--multiplier-degree", type=int, default=2)
    parser.add_argument("--atom-degree", type=int, default=2)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--compare-old", action="store_true")
    args = parser.parse_args()
    assert Path("/Volumes/data").is_mount()
    assert args.multiplier_degree >= 0 and args.atom_degree >= 1
    assert not args.cache_dir.exists(), "refusing to overwrite a matrix cache"
    I, J = parse_support(args.I), parse_support(args.J)
    Aeq, metadata = build(I, J, args.multiplier_degree, args.atom_degree, args.compare_old)
    args.cache_dir.mkdir(parents=True)
    matrix_path = args.cache_dir / "Aeq.int64.npy"
    np.save(matrix_path, Aeq, allow_pickle=False)
    metadata.update({
        "format": "positive-certificate-integer-matrix-cache-v1",
        "status": "matrix_built_no_lp_run",
        "source_sha256": {"builder": sha(Path(__file__)), "frozen_positive_search": sha(FROZEN)},
        "cache_file": matrix_path.name, "cache_file_sha256": sha(matrix_path),
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "sympy": S.__version__},
        "claim_boundary": "Exact integer LP-constraint matrix construction only. No LP was run and no certificate/nonexistence conclusion follows.",
    })
    (args.cache_dir / "matrix-metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: metadata[key] for key in ("matrix_shape", "free_multiplier_columns", "positive_atom_columns", "integer_matrix_sha256", "sparse_build_seconds", "old_matrix_comparison")}, sort_keys=True))


if __name__ == "__main__":
    main()
