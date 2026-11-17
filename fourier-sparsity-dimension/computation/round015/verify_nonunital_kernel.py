#!/usr/bin/env python3
"""Independent all-subspace verification of sigma(wB+zE)=kappa(B), n<=3."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time


def subspaces(n):
    for d in range(n+1):
        for piv in itertools.combinations(range(n), d):
            free = [(i, j) for i, p in enumerate(piv) for j in range(p+1, n) if j not in piv]
            for mask in range(1 << len(free)):
                b = [1 << p for p in piv]
                for k, (i, j) in enumerate(free):
                    if mask >> k & 1:
                        b[i] |= 1 << j
                yield b


def span(b):
    out = [0]
    for g in b:
        out += [x ^ g for x in out]
    return out


def capacity(cw, full):
    dp = {0: 0}
    for s in sorted(cw, key=int.bit_count):
        if s:
            dp[s] = max(sum(g & ~a == 0 for g in cw).bit_length() + dp[s ^ a]
                        for a in cw if a & (s & -s) and a & ~s == 0)
    return dp[full]


def main():
    p = argparse.ArgumentParser(); p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(); began = time.monotonic()
    rows = []; total_spaces = total_algebras = 0
    for n in range(1, 4):
        full = (1 << n)-1
        unital = [[full, *b] for b in subspaces(n-1)]
        for b in unital:
            bcw = span(b); bset = set(bcw); kb = capacity(bcw, full)
            for e in unital:
                eset = set(span(e))
                if not eset <= bset or not all((u & v) in eset for u in e for v in e):
                    continue
                if not all((u & v) in bset for u in e for v in b):
                    continue
                # Independent embedding on exactly 4*n points.
                nb = [sum(((g >> x) & 1) << (x+n*w+2*n*z)
                          for x in range(n) for w in range(2) for z in range(2) if w) for g in b]
                nb += [sum(((g >> x) & 1) << (x+n*w+2*n*z)
                           for x in range(n) for w in range(2) for z in range(2) if z) for g in e]
                nw = span(nb); nset = set(nw); assert len(nset) == 1 << len(nb)
                best = -1; checked = closed = 0; witness = None
                for fcoef in subspaces(len(nb)):
                    checked += 1
                    fw = {nw[u] for u in span(fcoef)}
                    fb = [nw[u] for u in fcoef]
                    if not all((u & v) in fw for u in fb for v in fb):
                        continue
                    closed += 1
                    # The maximal F-module in N. Closure follows from FF subset F;
                    # the check below also explicitly verifies it on all words.
                    dw = {d for d in nw if all((f & d) in nset for f in fb)}
                    assert fw <= dw and all((f & d) in dw for f in fb for d in dw)
                    assert len(dw) & (len(dw)-1) == 0
                    score = len(fcoef) + len(dw).bit_length()-1
                    if score > best:
                        best = score; witness = {'F_basis': fb, 'dim_D': len(dw).bit_length()-1}
                assert best == kb
                total_spaces += checked; total_algebras += closed
                rows.append(dict(points=n, B_basis=b, E_basis=e, kappa_B=kb,
                                 sigma_N=best, N_dimension=len(nb),
                                 subspaces_checked=checked, algebras_checked=closed,
                                 witness=witness))
    out = dict(status='passed', cases=len(rows), all_subspaces_checked=total_spaces,
               all_nonunital_algebras_checked=total_algebras, rows=rows,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               python=sys.version, elapsed_seconds=time.monotonic()-began,
               scope='All B,E pairs on n<=3 points satisfying 1 in E subset B and EB subset B; all linear subspaces F of each N, with maximal modules D. Finite check independent of the general proof.')
    args.out.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: out[k] for k in ['status','cases','all_subspaces_checked','all_nonunital_algebras_checked','elapsed_seconds']}))


if __name__ == '__main__':
    main()
