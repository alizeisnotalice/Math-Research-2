"""python3 -B check_components.py. Integer checks for matrix circuits and full spectral cuts."""
from fractions import Fraction
from itertools import combinations, product
from collections import Counter
import json
import platform


def vector_rank(columns):
    piv = {}
    for x in columns:
        while x and x.bit_length() in piv:
            x ^= piv[x.bit_length()]
        if x:
            piv[x.bit_length()] = x
    return len(piv)


def wedge(x,y,n):
    return sum((((x>>i)&1)*y ^ ((y>>i)&1)*x) << (i*n) for i in range(n))


def plane_catalog(n):
    result = {}
    for a,b in combinations(range(1,2**n),2):
        result[wedge(a,b,n)] = (a,b)
    return result


def components(forms):
    basis = []
    pivots = {}
    for f in forms:
        v, c = f, 0
        while v and v.bit_length() in pivots:
            w,u = pivots[v.bit_length()]
            v,c = v^w,c^u
        if v:
            c ^= 1<<len(basis)
            pivots[v.bit_length()] = v,c
            basis.append(f)
    groups = [{i} for i in range(len(basis))]
    for f in forms:
        v,c = f,0
        while v:
            w,u = pivots[v.bit_length()]
            v,c = v^w,c^u
        merged=set()
        for group in groups:
            if any(c>>i&1 for i in group):
                merged.update(group)
        groups = [g for g in groups if not (g & merged)] + [merged]
    return [[basis[i] for i in group] for group in groups]


def component_test():
    catalog = plane_catalog(4)
    stats = Counter()
    for size in range(1,5):
        for forms in combinations(catalog,size):
            for group in components(forms):
                b = vector_rank(x for a in group for x in catalog[a])
                assert b <= len(group)+1
                stats[(len(group),b)] += 1
    return {str(key):value for key,value in sorted(stats.items())}


def subspaces(n):
    for size in range(n+1):
        for pivots in combinations(range(n),size):
            free = [(i,column) for i,p in enumerate(pivots) for column in range(p)
                    if column not in pivots]
            for bits in range(1<<len(free)):
                rows = [1<<p for p in pivots]
                for j,(i,c) in enumerate(free):
                    rows[i] |= ((bits>>j)&1)<<c
                yield rows


def membership(x, rows):
    for r in reversed(rows):
        if x & (1 << (r.bit_length()-1)):
            x ^= r
    return x==0


def evaluate_seed(k):
    n = k+1
    pairs = [(1<<i,1<<(i+1)) for i in range(k)]
    pairs += [(1<<(i+1),(1<<i)^(1<<(i+2))) for i in range(k-1)]
    N = 1<<(len(pairs)-1).bit_length()
    m = N.bit_length()-1
    pairs += [pairs[0]]*(N-len(pairs))
    forms = [wedge(x,y,n) for x,y in pairs]
    groups = components(forms)
    assert [len(g) for g in groups] == [k]
    assert vector_rank([x for xy in pairs for x in xy]) == k+1
    hole = forms[0]^forms[1]^forms[2]
    hole_rank = vector_rank((hole>>(n*i))&((1<<n)-1) for i in range(n))
    if k>=5:
        hole = forms[0]^forms[2]^forms[4]
        hole_rank = vector_rank((hole>>(n*i))&((1<<n)-1) for i in range(n))
    assert hole_rank>=4
    old = max(Fraction((N+m-q-j+1)*N,N*N-q*(1<<j))
              for q in range(N+1) for j in range(m+1)
              if q<N or j<m)
    worst=Fraction(0)
    cuts=0
    for rows in subspaces(n):
        hits=sorted([1+sum(membership(z,rows) for z in (x,y,x^y)) for x,y in pairs],
                    reverse=True)
        sums=[0]
        for h in hits:
            sums.append(sums[-1]+h)
        for j in range(m+1):
            for q,inside in enumerate(sums):
                codim=N+m+n-q-j-len(rows)
                if not codim:
                    continue
                outside=4*N*N-(1<<j)*inside
                numerator=(codim+1)*2*N
                assert numerator*old.denominator<=outside*old.numerator
                if numerator*worst.denominator>outside*worst.numerator:
                    worst=Fraction(numerator,outside)
                cuts+=1
    return dict(path_edges=k, address_slots=N, joint_dimension=n, affine_hole_rank=hole_rank,
                old_L=str(old), new_L=str(worst), reduced_exact_cut_comparisons=cuts)


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),seed=None,tolerance=0,
                          matrix_checks=component_test(),
                          seed_families=[evaluate_seed(k) for k in (4,5,6)]),indent=2))
