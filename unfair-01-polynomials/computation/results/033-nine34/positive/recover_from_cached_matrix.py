#!/usr/bin/env python3
"""LP plus exact active-support recovery using ds_optimize's immutable integer cache."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from itertools import combinations_with_replacement, combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog
import sympy as S

ROOT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
HERE = Path(__file__).resolve().parent
CACHE = ROOT / 'agents/ds_optimize/round033-positive-fast/final-8x9-d2-a2-matrixcache'
TOPIC = Path('/Volumes/data/CodexProjects/T2/unfair-01-polynomials')
I = [0, 1, 2, 4, 6, 9, 11, 16]
J = [0, 1, 3, 4, 6, 8, 11, 13, 18]
NLEFT, NATOM = 2310, 338


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def pair(x: S.Rational) -> list[int]:
    x = S.Rational(x)
    return [int(x.p), int(x.q)]


def layouts() -> tuple[list[tuple[str, ...]], list[tuple[int, tuple[str, ...]]], list[tuple[tuple[str, ...], tuple[str, ...]]]]:
    names = ('a1', 'a2', 'a4', 'a6', 'a9', 'a11', 'b1', 'b3', 'b4', 'b6', 'b8', 'b11', 'b13')
    monos = [()] + [(x,) for x in names] + list(combinations_with_replacement(names, 2))
    diagonals = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 20, 22, 24, 27, 29)
    products = [(k, mono) for k in diagonals for mono in monos]
    atoms = []
    for r in (1, 2):
        for picks in combinations(range(len(names)), r):
            for signs in product((0, 1), repeat=r):
                atoms.append((tuple(names[t] for t, sign in zip(picks, signs) if sign == 0),
                              tuple(names[t] for t, sign in zip(picks, signs) if sign == 1)))
    assert len(products) == NLEFT and len(atoms) == NATOM
    return list(names), products, atoms


def main() -> None:
    global HERE,CACHE,TOPIC
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);parser.add_argument('--cache',type=Path,required=True);parser.add_argument('--topic',type=Path,required=True);args=parser.parse_args()
    HERE,CACHE,TOPIC=args.out,args.cache,args.topic
    assert Path('/Volumes/data').is_mount()
    HERE.mkdir(parents=True,exist_ok=True)
    capture, candidate, verify = HERE/'cached-active-basis-capture.json', HERE/'cached-recovered-candidate.json', HERE/'cached-independent-verification.json'
    assert not any(p.exists() for p in (capture, candidate, verify)), 'refusing overwrite'
    meta = json.loads((CACHE/'matrix-metadata.json').read_text())
    assert meta['I'] == I and meta['J'] == J and meta['integer_matrix_sha256'] == 'a07064e462c3bbb2586b815893441278c14a35a3ab0cf1c6e6c21d880544c15d'
    assert sha(CACHE/'Aeq.int64.npy') == meta['cache_file_sha256']
    A = np.load(CACHE/'Aeq.int64.npy', mmap_mode='r')
    assert A.shape == (2381, 4958) and A.dtype == np.dtype('<i8')
    b = np.zeros(A.shape[0]); b[-1] = 1
    answer = linprog(np.r_[np.ones(2*NLEFT), np.zeros(NATOM)], A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    assert answer.success, answer.message
    active = np.flatnonzero(answer.x > 1e-9)
    raw = {'format': 'round033-cached-active-basis-v1', 'matrix_cache': {'path': str(CACHE), 'metadata_sha256': sha(CACHE/'matrix-metadata.json'),
           'file_sha256': sha(CACHE/'Aeq.int64.npy'), 'integer_matrix_sha256': meta['integer_matrix_sha256'], 'shape': list(A.shape)},
           'solver': {'name': 'SciPy HiGHS linprog', 'message': answer.message, 'objective': float(answer.fun)},
           'active_threshold': 1e-9, 'active_columns': [{'index': int(i), 'float_hex': float(answer.x[i]).hex(), 'value': float(answer.x[i])} for i in active]}
    capture.write_text(json.dumps(raw, ensure_ascii=False, indent=2)+'\n')
    # Numeric QR only selects an independent row set; all final equalities below use exact integers.
    active_matrix = np.asarray(A[:, active], dtype=float)
    _, r, pivots = qr(active_matrix.T, pivoting=True, mode='economic')
    rank = int(np.sum(np.abs(np.diag(r)) > 1e-9))
    selected_rows = [int(x) for x in pivots[:rank]]
    raw['exact_recovery'] = {'active_count': int(len(active)), 'numeric_rank': rank, 'selected_independent_rows': selected_rows}
    if rank != len(active):
        raw['exact_recovery']['status'] = 'active_support_not_independent; no arbitrary rational point selected'
        capture.write_text(json.dumps(raw, ensure_ascii=False, indent=2)+'\n')
        print(json.dumps(raw['exact_recovery'], ensure_ascii=False)); return
    M = S.Matrix(np.asarray(A[selected_rows][:, active], dtype=np.int64).tolist())
    rhs = S.Matrix([int(b[row]) for row in selected_rows])
    exact = [S.Rational(x) for x in M.inv() * rhs]
    full = S.Matrix(np.asarray(A[:, active], dtype=np.int64).tolist()) * S.Matrix(exact)
    exact_ok = full == S.Matrix([int(x) for x in b]) and all(x >= 0 for x in exact)
    raw['exact_recovery'].update({'exact_nonnegative_full_matrix_check': exact_ok, 'min_value': str(min(exact))})
    if not exact_ok:
        raw['exact_recovery']['status'] = 'active-basis exact solve failed nonnegativity or full equality'
        capture.write_text(json.dumps(raw, ensure_ascii=False, indent=2)+'\n'); print(json.dumps(raw['exact_recovery'], ensure_ascii=False)); return
    names, products, atoms = layouts()
    z = [S.Integer(0)] * A.shape[1]
    for column, value in zip(active, exact): z[int(column)] = value
    p = [z[k] - z[NLEFT+k] for k in range(NLEFT)]
    q = z[2*NLEFT:]
    record = {'format': 'positive-certificate-v1', 'case': {'name': 'nine34-m16-line11-cache-active-basis', 'I': I, 'J': J},
      'strict_variables': names, 'unit_variables': ['a0', 'a16', 'b0', 'b18'],
      'basis': {'multiplier_total_degree_at_most': 2, 'positive_atom_total_degree_at_most': 2, 'atom_form': 'product(positive) * product(1-minus)'},
      'certificate': {'multipliers': [{'diagonal': k, 'monomial': list(mono), 'coefficient': pair(v)} for v,(k,mono) in zip(p,products) if v],
       'positive_atoms': [{'positive': list(pos), 'one_minus': list(neg), 'coefficient': pair(v)} for v,(pos,neg) in zip(q,atoms) if v], 'exact_recovery_check': True},
      'status': 'certificate_found'}
    candidate.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
    run = subprocess.run([sys.executable, str(TOPIC/'computation/verify_positive_certificate.py'), str(candidate)], text=True, capture_output=True)
    verify.write_text(json.dumps({'returncode': run.returncode, 'stdout':run.stdout, 'stderr':run.stderr}, ensure_ascii=False, indent=2)+'\n')
    raw['exact_recovery']['independent_fraction_verifier_returncode'] = run.returncode
    raw['exact_recovery']['status'] = 'certificate_verified' if run.returncode == 0 else 'candidate_rejected_by_independent_verifier'
    capture.write_text(json.dumps(raw, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(raw['exact_recovery'], ensure_ascii=False))


if __name__ == '__main__':
    main()
