#!/usr/bin/env python3
"""Independent capacity checks of all 45 base completions by algebra enumeration."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
from build_extension_seed import basis, span


def subspaces(n):
    for d in range(n+1):
        for piv in itertools.combinations(range(n), d):
            free = [(i, j) for i, p in enumerate(piv) for j in range(p+1, n) if j not in piv]
            for m in range(1 << len(free)):
                b = [1 << p for p in piv]
                for k, (i, j) in enumerate(free):
                    if m >> k & 1:
                        b[i] |= 1 << j
                yield b


def check(item):
    index, original = item
    adapted = [65535]
    for v in original:
        if len(basis([*adapted, v])) > len(adapted):
            adapted.append(v)
    assert len(adapted) == 7
    code = span(adapted); cset = set(code)
    tested = closed = 0; best = 0; witness = None
    for qb in subspaces(6):
        tested += 1
        # Parameter bit 0 is the full-domain unit, bits 1..6 the quotient.
        fb = [65535, *(code[v << 1] for v in qb)]
        fset = set(span(fb))
        if not all((a & b) in fset for a in fb for b in fb):
            continue
        closed += 1
        core = [d for d in code if all((f & d) in cset for f in fb)]
        cdim = len(core).bit_length()-1
        assert len(core) == 1 << cdim
        assert all((f & d) in set(core) for f in fb for d in core)
        score = len(fb) + cdim
        if score > best:
            best = score; witness = {'algebra_basis': fb, 'module_basis': list(basis(core))}
    assert tested == 2825 and best == 8
    return dict(index=index, unital_subspaces=tested, algebras=closed,
                kappa=best, witness=witness)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--input', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True); p.add_argument('--workers', type=int, default=12)
    a = p.parse_args(); start = time.monotonic()
    source = json.loads(a.input.read_text()); completions = source['maximum_base_extension_bases']
    assert len(completions) == 45
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(check, enumerate(completions)))
    out = dict(status='passed', cases=len(rows), workers=a.workers,
               total_unital_subspaces=sum(r['unital_subspaces'] for r in rows),
               total_algebras=sum(r['algebras'] for r in rows), rows=rows,
               input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               utility_sha256=hashlib.sha256(Path(__file__).with_name('build_extension_seed.py').read_bytes()).hexdigest(),
               python=sys.version, elapsed_seconds=time.monotonic()-start,
               scope='Independent enumeration of all unital subalgebras and their maximal modules, not the support DP used by the seed builder.')
    a.out.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','cases','workers','total_unital_subspaces','total_algebras','elapsed_seconds']}))


if __name__ == '__main__':
    main()
