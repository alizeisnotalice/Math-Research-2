#!/usr/bin/env python3
"""Recompute strict instances using set partitions, rank checks and the R13 DP."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from triangle_graph_scan import basis, words


def partitions(n):
    def visit(cells, i):
        if i == n:
            yield cells
        else:
            for j in range(len(cells)):
                yield from visit([c | (1 << i) if k == j else c for k, c in enumerate(cells)], i+1)
            yield from visit([*cells, 1 << i], i+1)
    yield from visit([], 0)


def kap_by_set_partitions(b, n):
    code = set(words(b)); d = len(b)
    best = -1; witness = None; tried = 0
    for p in partitions(n):
        tried += 1
        if all(c in code for c in p):
            value = sum(d-len(basis(g & ~c for g in b))+1 for c in p)
            if value > best:
                best = value; witness = p
    return {'kappa': best, 'partition': witness, 'all_set_partitions': tried}


def main():
    p = argparse.ArgumentParser()
    for name in ['inputs', 'binary', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(); began = time.monotonic(); raw = json.loads(a.inputs.read_text())
    rows = []
    for name, ex in raw['examples'].items():
        cb, bb, e = ex['C_basis'], ex['B_basis'], ex['E_cells']
        cset, bset = set(words(cb)), set(words(bb))
        assert set(words(e)) <= bset <= cset
        assert all((u & v) in bset for u in e for v in bb)
        kc = kap_by_set_partitions(cb, 7); kb = kap_by_set_partitions(bb, 7)
        assert kc['kappa'] == ex['kappa_C'] and kb['kappa'] == ex['kappa_B']
        def lift(g, part):
            return sum(((g >> (i % 28 % 7)) & 1) << i for i in range(64)
                       if part == 'C' or (part == 'B' and ((i % 28)//7) & 1)
                       or (part == 'E' and ((i % 28)//14) & 1))
        tb = basis([*(lift(g, 'C') for g in cb), *(lift(g, 'B') for g in bb), *(lift(g, 'E') for g in e)])
        assert tb == basis(ex['T_lift_basis'])
        d = len(tb)
        run = subprocess.run([str(a.binary.resolve()), str(d), *map(str, tb)], capture_output=True, text=True, timeout=30, check=True)
        result = json.loads(run.stdout)
        assert result['D'] == d and result['kappa'] == ex['kappa_T']
        # Check both the historical partition and the fresh one.
        for part in [ex['T_optimal_partition'], result['partition']]:
            union = score = 0
            for s in part:
                cell = int(s, 16)
                assert cell and not cell & union and len(basis([*tb, cell])) == d
                union |= cell; score += d-len(basis(g & ~cell for g in tb))+1
            assert union == (1 << 64)-1 and score == ex['kappa_T']
        rows.append(dict(name=name, C=kc, B=kb, T=result, exit_code=run.returncode,
                         gap=kc['kappa']+kb['kappa']-result['kappa']))
    out = dict(status='passed', rows=rows, python=sys.version,
               elapsed_seconds=time.monotonic()-began,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               utility_sha256=hashlib.sha256(Path(__file__).with_name('triangle_graph_scan.py').read_bytes()).hexdigest(),
               input_sha256=hashlib.sha256(a.inputs.read_bytes()).hexdigest(),
               binary_sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
               scope='Fresh local computations: all 877 set partitions of each C/B, full integer DP of each T, and both historical/fresh T partition rank certificates. No bentness assertion for these abstract function codes.')
    a.out.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status':'passed','examples':len(rows),'gaps':[r['gap'] for r in rows],'elapsed_seconds':out['elapsed_seconds']}))


if __name__ == '__main__':
    main()
