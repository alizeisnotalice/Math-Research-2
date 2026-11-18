"""Exact correlations for all three-variable quadratic pairs, modulo constant signs."""
import json
import platform
from collections import Counter
from itertools import combinations, product


def bit(x, k):
    return (x >> k) & 1


def walsh(q, n):
    return [sum((-1)**(bit(q,x) ^ ((x&z).bit_count()%2)) for x in range(1<<n))
            for z in range(1<<n)]


def rank(q, n):
    roots = sum(all(not (bit(q,0)^bit(q,x)^bit(q,y)^bit(q,x^y)) for y in range(1<<n))
                for x in range(1<<n))
    return n - roots.bit_length() + 1


def audit(n, include_linear):
    N = 1 << n
    monomials = [1 << j for j in range(n)] if include_linear else []
    monomials += [(1 << j) | (1 << k) for j,k in combinations(range(n),2)]
    laws = {0}
    for m in monomials:
        term = sum(1 << x for x in range(N) if (x&m)==m)
        laws |= {q ^ term for q in laws}
    spectra = {q: walsh(q,n) for q in laws}
    roots = {q: rank(q,n) for q in laws}
    disjoint = 0
    types = Counter()
    for f,g in product(laws,repeat=2):
        P,Q=spectra[f],spectra[g]
        I=[z for z in range(N) if P[z] and Q[z]]
        if not I:
            disjoint += 1
            continue
        c=max(abs(sum((1 if P[z]*Q[z]>0 else -1)*(-1)**((z&t).bit_count()%2) for z in I))
              for t in range(N))
        a,b=roots[f],roots[g]
        rho=roots[f^g]
        assert c*c*(1<<rho) == 1<<(a+b)
        lam=2*(len(I).bit_length()-c.bit_length())
        i=len(I).bit_length()-1
        assert rho==(a-i)+(b-i)+lam
        types[(a-i,b-i,lam)] += 1
    return {"ordered_pairs":len(laws)**2, "disjoint":disjoint, "nonempty_intersection":sum(types.values()),
            "cases_with_unequal_spectra":sum(count for (u,v,_),count in types.items() if u or v),
            "private_codimensions_and_common_phase_rank":{
                ",".join(map(str,z)):t for z,t in sorted(types.items())}}

print(json.dumps({"python":platform.python_version(), "all_three_variable_pairs":audit(3,True),
                  "four_variable_homogeneous_pairs":audit(4,False)}, indent=2))
