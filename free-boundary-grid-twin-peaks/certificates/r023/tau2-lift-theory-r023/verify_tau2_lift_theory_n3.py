#!/usr/bin/env python3
"""Exact r023 checks: interleaving tails/moments and n3 first-merge fibres."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path

N = 9
HERE = Path(__file__).resolve().parent


def xy(v): return v % 3, v // 3


ADJ = [set() for _ in range(N)]
for v in range(N):
    x, y = xy(v)
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        X, Y = x + dx, y + dy
        if 0 <= X < 3 and 0 <= Y < 3:
            ADJ[v].add(Y * 3 + X)


def births(word):
    seen, ans = set(), []
    for r, v in enumerate(word, 1):
        if not (seen & ADJ[v]): ans.append((r, v))
        seen.add(v)
    return ans


def components(vertices):
    left, ans = set(vertices), []
    while left:
        todo, c = [left.pop()], set()
        while todo:
            v = todo.pop(); c.add(v)
            new = ADJ[v] & left; left -= new; todo.extend(new)
        ans.append(frozenset(c))
    return ans


def fibre(word):
    tau = births(word)[1][0]
    for r in range(tau + 1, N + 1):
        old, w = set(word[:r - 1]), word[r - 1]
        if sum(bool(ADJ[w] & c) for c in components(old)) < 2: continue
        S = frozenset(word[:r]); A, B = sorted(components(set(S)-{w}), key=lambda c:min(c))
        seqA = tuple(v for v in word[:r-1] if v in A)
        seqB = tuple(v for v in word[:r-1] if v in B)
        return (tuple(sorted(S)), w, seqA, seqB, tuple(word[r:])), len(A), len(B)
    raise AssertionError


def L(word):
    y = word[1]
    return next(r for r, v in enumerate(word, 1) if r >= 3 and v in ADJ[y])-2


def e2(word):
    x, y = xy(word[0]), xy(word[1])
    return (x[0]-y[0])**2+(x[1]-y[1])**2


def ab_samples(a, b):
    q, ans = a+b, []
    for slots in combinations(range(q-2), a-1):
        letters = ['B']*(q-2)
        for i in slots: letters[i] = 'A'
        ans.append(q-1 if 'B' not in letters else letters.index('B')+1)
    return ans


def second(a,b):
    return Fraction(1) + Fraction(3*(a-1),b) + Fraction(2*(a-1)*(a-2),b*(b+1))


def main():
    small = 0
    for a in range(1, 13):
        for b in range(1, 13):
            xs = ab_samples(a,b)
            mean, sq = Fraction(sum(xs),len(xs)), Fraction(sum(x*x for x in xs),len(xs))
            assert mean == Fraction(a+b-1,b)
            assert sq == second(a,b)
            assert sq <= 2*mean*mean
            for ell in range(1,a+2):
                observed = Fraction(sum(x >= ell for x in xs),len(xs))
                expected = (Fraction(0) if ell > a else
                            Fraction(__import__('math').comb(a-1,ell-1), __import__('math').comb(a+b-2,ell-1)))
                assert observed == expected
            small += 1

    by_fibre = defaultdict(list)
    for word in permutations(range(N)):
        b = births(word)
        if len(b) == 2 and b[1][0] == 2:
            key,a,bb = fibre(word)
            by_fibre[key].append((L(word), a, bb, e2(word)))
    sumL2 = sumM = Fraction(); sumM2 = Fraction(); sumLD = 0; sumMD = Fraction()
    for values in by_fibre.values():
        a, bb = values[0][1:3]
        assert all((u,v)==(a,bb) for _,u,v,_ in values)
        q = a+bb
        m = Fraction(q*(q-1),2*a*bb)
        s = (second(a,bb)+second(bb,a))/2
        assert Fraction(sum(z*z for z,_,_,_ in values),len(values)) == s
        assert s <= 4*m*m
        assert Fraction(sum(z for z,_,_,_ in values),len(values)) == m
        # Root distance is fixed by the fibre data, hence so is H.
        assert len({d for *_,d in values}) == 1
        d = values[0][3]
        sumL2 += sum(z*z for z,_,_,_ in values)
        sumM += len(values)*m; sumM2 += len(values)*m*m
        sumLD += sum(z*d for z,_,_,d in values)
        sumMD += len(values)*m*d
    assert sumM == 143808
    assert sumLD == 633600 and sumMD == 633600
    assert sumL2 == 358080 and sumL2 <= 4*sumM2
    out = {"n":3, "method":"all tau=2 words from all 9! permutations; exact Fraction fibre moments",
           "small_ab_pairs_checked":small,"complete_fibres":len(by_fibre),
           "sum_L2":str(sumL2),"sum_m":str(sumM),"sum_m2":str(sumM2),
           "sum_LdE2":sumLD,"sum_mdE2":str(sumMD),
           "checks":{"all_tail_formulas":True,"all_direction_moment_formulas":True,
                     "all_fibre_RB_moment_formula":True,"global_L2_le_4_m2":True,
                     "RB_euclidean_total":True}}
    (HERE/'tau2-lift-theory-n3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
