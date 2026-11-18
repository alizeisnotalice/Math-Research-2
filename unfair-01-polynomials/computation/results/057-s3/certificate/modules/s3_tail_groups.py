#!/usr/bin/env python3
"""Stream exact S={3} Newton tail groups and optionally intersect Cayley halfplanes."""

from __future__ import annotations

import argparse
import hashlib
from math import comb, isqrt
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from s3_entry_frozen import (
    BOUND, INITIAL, ONE, ZERO, ci, cores, count_translated_row, integer_interval,
    newton_next, prepare_translated_row, row_step, strict_house, strict_integer_ball,
    newton_without_current, add, tail_box,
)


def comb_or_zero(n, k):
    return comb(n,k) if 0 <= k <= n else 0


# C[j][k] = [y^k] 2^j (1+y)^j (1-y)^(8-j).
CAYLEY = tuple(tuple(2**j * sum((-1)**(k-i)*comb_or_zero(j,i)*comb_or_zero(8-j,k-i)
                             for i in range(max(0,k-(8-j)),min(j,k)+1))
                      for k in range(9)) for j in range(9))


def cayley_inequalities(L, T, D, b, c):
    """Return (D coefficient, E coefficient, upper RHS) for lower halfspace."""
    co = (2, 2*L, 2*L*L, T, 2*b, 2*c, 0, 0, 2)
    result = []
    for reverse in (False, True):
        mat = CAYLEY[::-1] if reverse else CAYLEY
        for k in range(9):
            base = sum(x * column[k] for x, column in zip(co, mat, strict=True))
            minimum = isqrt(mat[3][k]**2*D) + 1
            result.append((2*mat[6][k], 2*mat[7][k], minimum-base))
    return result


def grouped_intervals(core):
    """Yield b,c and exact constraints on d,e from low prefix and Newton 1--8."""
    L, M, T, N, D = core
    p = [ONE,ci(L),ci(M),(T,1),ZERO,ZERO,ZERO,ZERO,ONE]
    s = [ci(8)]
    states = INITIAL
    for k in range(1, 4):
        s.append(newton_next(p,s,k,D))
        states = row_step(states,p,k,D)
    base4 = newton_without_current(p,s,4,D)
    for b in strict_integer_ball(base4,4,D):
        p[4] = ci(b)
        s4 = add(base4,ci(-4*b))
        live = states
        for k in range(4,8):
            live = row_step(live,p,k,D)
            if not live:
                break
        if not live:
            continue
        prepared = prepare_translated_row(live,p,8,D)
        base5 = newton_without_current(p,[*s,s4],5,D)
        for c in strict_integer_ball(base5,5,D):
            if not count_translated_row(prepared,c,D):
                continue
            p[5] = ci(c)
            s5 = add(base5,ci(-5*c))
            powers = [*s,s4,s5]
            for k in range(6,9):
                powers.append(newton_next(p,powers,k,D))
            elo, ehi = -9, 9
            dlo, dhi = integer_interval(powers[6],BOUND[6],-6,D)
            dlo, dhi = max(-41,dlo), min(41,dhi)
            zlo, zhi = integer_interval(powers[7],BOUND[7],-7,D)
            if L:
                elo8, ehi8 = integer_interval(powers[8],BOUND[8],-8*abs(L),D)
                if L > 0:
                    elo8, ehi8 = -ehi8,-elo8
                elo, ehi = max(elo,elo8), min(ehi,ehi8)
            elif not strict_house(powers[8],*BOUND[8],D):
                continue
            if dlo > dhi or elo > ehi or zlo > zhi:
                continue
            result = tail_box(p, [*s,s4,s5], D)
            if result is not None:
                yield b,c,dlo,dhi,elo,ehi,zlo,zhi,result[0]
        p[5] = ZERO


def main():
    if sys.flags.optimize:
        raise SystemExit('Python assertions must be enabled')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--samples-per-L',type=int,help='deterministic, evenly spaced pilot within each L')
    ap.add_argument('--core-offset',type=int,default=0)
    ap.add_argument('--core-limit',type=int)
    args = ap.parse_args()
    assert args.core_limit is None or args.samples_per_L is None
    all_cores = cores()
    if args.samples_per_L is not None:
        assert args.samples_per_L > 0
        use_cores = []
        for L in range(-4,5):
            block = [x for x in all_cores if x[0]==L]
            use_cores.extend(block[(2*i+1)*len(block)//(2*args.samples_per_L)] for i in range(args.samples_per_L))
    elif args.core_limit is not None:
        use_cores = all_cores[args.core_offset:args.core_offset+args.core_limit]
    else:
        use_cores = all_cores
    hash_groups = hashlib.sha256()
    n_groups = n_mothers = 0
    with args.out.open('xb') as out:
        for core in use_cores:
            L, _, T, N, D = core
            for fields in grouped_intervals(core):
                n_groups += 1
                n_mothers += fields[-1]
                line = (' '.join(map(str, (L,T,N,D,*fields)))+'\n').encode()
                hash_groups.update(line)
                out.write(line)
    print(f'cores={len(use_cores)} groups={n_groups} mothers={n_mothers} sha256={hash_groups.hexdigest()}')


if __name__ == '__main__':
    main()
