"""Exact checks of Boolean clipping, quotient copies, and cyclic spanning repair."""

import json
import random
import sys
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import comb

import numpy as np
import scipy
from scipy.optimize import linprog

SEED = 860037


@cache
def span(vectors):
    space = {0}
    for vector in vectors:
        if vector not in space:
            space |= {point ^ vector for point in space}
    return frozenset(space)


@cache
def all_subspaces(space):
    known = {frozenset([0])}
    pending = list(known)
    for part in pending:
        for v in space-part:
            enlarged = part | frozenset(p^v for p in part)
            if enlarged not in known:
                known.add(enlarged)
                pending.append(enlarged)
    return pending


def walsh_magnitudes(table):
    result = []
    for frequency in range(8):
        total = sum((-1 if ((index & frequency).bit_count() % 2) else 1)*value
                    for index,value in enumerate(table))
        if total:
            result.append((8+frequency,abs(total)))
    return tuple(result)


def quotient(parent,H):
    return tuple(sorted((representative,magnitude) for vector,magnitude in parent
                        if (representative := min(vector^h for h in H))))


def rank(vectors):
    return len(span(tuple(sorted(vectors)))).bit_length()-1


def polytope_distribution(candidates, capacities, optimize_scale=False):
    count = len(capacities)
    incidence = [[int(e in selected) for selected in candidates] for e in range(count)]
    if optimize_scale:
        A = [row + [-float(cap)] for row,cap in zip(incidence,capacities)]
        B = [0.0]*count
        M = [1.0]*len(candidates)+[0.0]
        objective = [0.0]*len(candidates)+[1.0]
    else:
        A = incidence
        B = [float(cap) for cap in capacities]
        M = [1.0]*len(candidates)
        objective = [0.0]*len(candidates)
    result = linprog(objective,A_ub=A,b_ub=B,A_eq=[M],b_eq=[1.0],
                     bounds=(0,None),method="highs",options={"threads":1})
    assert result.success
    weights = [
        (selected,Fraction(float(weight)).limit_denominator(100_000))
        for selected,weight in zip(candidates,result.x) if weight>1e-9
    ]
    total = sum(weight for _,weight in weights)
    return [(selected,weight/total) for selected,weight in weights]


def coalition_repair(active, columns, raw_weights, prices):
    whole_rank = rank(columns)
    bad = [e for e in active if rank([columns[j] for j in active if j!=e])<whole_rank]
    generator = {e:span(tuple(columns[j] for j in active if j!=e)) for e in bad}
    covers = {f:frozenset(e for e in bad if columns[f] not in generator[e])
              for f in range(len(columns)) if f not in active}
    assert all(sum(raw_weights[f] for f in covers if e in covers[f]) >= Fraction(1,2)
               for e in bad)
    unresolved = set(bad)
    repairs = set()
    while unresolved:
        candidates = [f for f in covers if covers[f]&unresolved]
        best = min(candidates,key=lambda f:Fraction(prices[f],len(unresolved&covers[f])))
        repairs.add(best)
        unresolved -= covers[best]
    selected = set(active)|repairs
    assert all(rank([columns[j] for j in selected if j!=e])==whole_rank for e in selected)
    return sum(prices[e] for e in selected),len(bad)


def check(instance, rng):
    columns = tuple(point for point,_ in instance)
    weights = tuple(Fraction(value,8) for _,value in instance)
    V = span(columns)
    r = len(V).bit_length()-1
    proper = [U for U in all_subspaces(V) if U!=V]
    L = max((r+2-len(U).bit_length())/sum(w for col,w in zip(columns,weights) if col not in U)
            for U in proper)
    capacities=[L*w for w in weights]
    for U in proper:
        q = r+1-len(U).bit_length()
        assert sum(min(b,1) for p,b in zip(columns,capacities) if p not in U)>=q+1
    for test in range(16):
        exposed = [w for e,w in zip(columns,weights) if (test&e).bit_count()%2]
        if exposed:
            assert len(exposed)>1 and sum(exposed)>=1
            assert sum(exposed) - max(exposed)>=Fraction(1,2)

    ids = list(range(len(columns)))
    bases = [subset for subset in combinations(ids,r) if rank([columns[i] for i in subset])==r]
    all_cyclic = []
    for size in range(r+1,len(ids)+1):
        for subset in combinations(ids,size):
            if all(rank([columns[i] for i in subset if i != removed])==r for removed in subset):
                all_cyclic.append(subset)
    base_distribution = polytope_distribution(bases,capacities)
    probabilities = [sum(weight for basis,weight in base_distribution if i in basis) for i in ids]
    assert all(p<=b for p,b in zip(probabilities,capacities))
    cyclic_distribution = polytope_distribution(all_cyclic,capacities,True)
    ratio = max(sum(weight for subset,weight in cyclic_distribution if i in subset)/capacities[i]
                for i in ids)
    assert ratio>=1
    if ratio>2:
        assert (ratio-2)**2 <= Fraction(9,2)

    prices=[rng.randrange(1,16) for _ in ids]
    expected_cost=Fraction(0)
    expected_h=Fraction(0)
    for basis_a,p in base_distribution:
        for basis_b,q in base_distribution:
            active = frozenset(basis_a)|frozenset(basis_b)
            cost,h=coalition_repair(active,columns,weights,prices)
            expected_cost += p*q*cost
            expected_h += p*q*h
    assert expected_h<=sum(p**2 for p in probabilities)<=L**2
    normalized_cost = sum(b*c for b,c in zip(capacities,prices))
    if expected_cost > 2*normalized_cost:
        assert (expected_cost-normalized_cost*2)**2 <= Fraction(9,2)*normalized_cost**2
    return L,ratio


def main():
    rng=random.Random(SEED)
    unique={walsh_magnitudes(table) for table in product((-1,1),repeat=8)}
    noncharacters={instance for instance in unique if len(instance)>1}
    quotient_instances=set()
    original_cuts=0
    for parent in noncharacters:
        V=span(tuple(col for col,_ in parent))
        for space in all_subspaces(V):
            projected=quotient(parent,space)
            if projected:
                quotient_instances.add(projected)
                original_cuts += 1

    best=Fraction(0)
    verified=0
    for instance in sorted(quotient_instances):
        L,ratio=check(instance,rng)
        best=max(best,ratio)
        verified+=1
    # All masks of k chosen coordinate sets omit a disjoint k-subset.
    setcover_cases=[]
    for k in (2,4,8):
        count=comb(2*k-1,k-1)
        assert count%2==1
        for chosen in combinations(range(2*k),k):
            available=set(range(2*k))-set(chosen)
            assert len(available)==k
        setcover_cases.append({"k":k,"cheap_basis_size":comb(2*k,k),
                               "distinct_covers_needed":k+1,
                               "abstract_multiplier_bound":str(Fraction(k+1,2))})
    print(json.dumps({"python":sys.version.split()[0],"scipy":scipy.__version__,
                      "seed":SEED,"distinct_noncharacter_weighted_spectra":len(noncharacters),
                      "actual_parent_and_quotient_cases":original_cuts,
                      "distinct_tested_quotient_instances":verified,
                      "max_rational_feasible_cyclic_ratio_from_lp":str(best),
                      "abstract_set_cover_cases":setcover_cases},indent=2))


if __name__=="__main__":
    main()
