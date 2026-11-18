#!/usr/bin/env python3
"""Find a real Boolean restriction where separated exterior spans meet in the quotient."""
from collections import Counter
from itertools import combinations
import json
import sys

import check_components as core


def parity(x):
    return bin(x).count('1') % 2


def spectrum(truth,n):
    return {u:a for u in range(1<<n)
            if (a:=sum(1-2*(((truth>>x)&1)^parity(x&u)) for x in range(1<<n)))}


def span(points):
    space={0}
    for x in points:
        space.update(x^y for y in list(space))
    return space


def project_coefficients(raw,ground,pivots,x):
    values=Counter()
    for f in ground:
        values[core.reduce_by(f,pivots)] += raw[f]*(1-2*parity(f&x))
    return {f:a for f,a in sorted(values.items()) if f and a}


def find(n, min_active_labels=1, n_deleted=(1,2)):
    size=1<<n
    for truth in range(1<<(size-1)):
        raw=spectrum(truth,n)
        if len(raw)<4:
            continue
        support=sum(1<<w for w in raw)
        for codim,H in core.spaces(n):
            if codim<2:
                continue
            exterior=support&~H
            for count in n_deleted:
                for deleted in combinations(core.points(exterior),count):
                    K=span([*core.points(H),*deleted])
                    if len(K)==len(set(core.points(H))) or len(K)==size:
                        continue
                    remain=exterior^sum(1<<u for u in deleted)
                    cs=core.compute_components(n,remain)
                    if len(cs)!=2:
                        continue
                    p=core.h_reducer(sum(1<<u for u in K))
                    groups=[list(core.points(c)) for c in cs]
                    original_spans=[span(u | size for u in group) for group in groups]
                    if original_spans[0]&original_spans[1]!={0}:
                        raise AssertionError('components not in direct sum')
                    for x in range(size):
                        active=[project_coefficients(raw,group,p,x) for group in groups]
                        if not all(active) or len(set(active[0]) | set(active[1])) < min_active_labels:
                            continue
                        overlap=set(active[0])&set(active[1])
                        witness=[q for q in sorted(overlap) if active[0][q]+active[1][q]]
                        if not witness:
                            continue
                        q=witness[0]
                        direct=project_coefficients(raw,list(raw),p,x)
                        if direct[q]!=active[0][q]+active[1][q]:
                            raise AssertionError('contributions do not add up in true slice')
                        leaf_points=sorted(y for y in range(size) if all(not parity((x^y)&u) for u in K))
                        direct_leaf=sum((1-2*((truth>>y)&1))*(1-2*parity(q&(x^y))) for y in leaf_points)
                        if direct_leaf*len(K)!=direct[q]:
                            raise AssertionError('independent restricted Walsh calculation failed')
                        return dict(n=n,truth_mask=truth,integer_walsh=raw,H=list(core.points(H)),removed=list(deleted),
                                    K=sorted(K),components=groups,projected_active_coefficients=active,
                                    common_nonzero_label=q,real_slice_coefficients=direct,
                                    real_slice_shift=x,real_slice_points=leaf_points,
                                    real_slice_values=[1-2*((truth>>y)&1) for y in leaf_points])
    return None


if __name__=='__main__':
    w=next((out for n in (2,3,4) if (out:=find(n)) is not None),None)
    if w is None:
        sys.exit('no quotient overlap found in tested range')
    stronger=next((out for n in (3,4) if (out:=find(n,3,(1,))) is not None),None)
    if stronger is None:
        sys.exit('no higher-rank quotient overlap in tested range')
    print(json.dumps(dict(status='pass',witness=w,stronger_witness=stronger),indent=2,sort_keys=True))
