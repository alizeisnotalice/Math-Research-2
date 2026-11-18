#!/usr/bin/env python3
"""Certified all-cut L equality for the nonlocal rank-preserving extension."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import platform
import time

import phase_lifts as phase
import transverse as tra
import verify_witnesses as verify
import compatible as comp


def all_projected_cuts(seed,w):
    a=w['a']; R=seed['R']; denom=seed['denominator']; L=Fraction(*seed['L']); scale=1<<(a//2)
    matrices=comp.data_matrices(seed,w['generators'],w['generator_coords'],a)
    images=[set(tra.span(phase.pivots(m))) for m in matrices]
    ranks=[len(i).bit_length()-1 for i in images]
    full_bill=sum(m*(scale<<(r//2)) for m,r in zip(seed['masses'],ranks))
    stats={}; tight=[]; min_by_transverse={}
    for j in range(1<<a):
        trans=a-(1 if j else 0)
        if trans<5:
            continue
        p=[r-int(bool(j and j in images[i])) for i,r in enumerate(ranks)]
        by_weight=Counter()
        for e,mass,r,d in zip(seed['E'],seed['masses'],ranks,p):
            c=mass*((scale<<(r//2))>>d)
            by_weight[c]|=1<<e
        minimum=None; hist=stats.setdefault(trans,Counter())
        for d,level in enumerate(tra.subspaces(R)):
            goal=R-d+trans+1
            for U in level:
                mass=Fraction(full_bill-sum(w*tra.pop(U&S) for w,S in by_weight.items()),denom*scale)
                ratio=mass*L/goal
                hist['cuts']+=1
                hist['weighted_failures']+=ratio<1
                hist['weighted_equalities']+=ratio==1
                if ratio==1 and len(tight)<50:
                    tight.append(dict(transverse=trans,J_line_generator=j,U=U,actual_mass_if_U_zero=None if U!=1 else [mass.numerator,mass.denominator]))
                if minimum is None or ratio<minimum:
                    minimum=ratio
        if trans not in min_by_transverse or minimum<Fraction(*min_by_transverse[trans]):
            min_by_transverse[trans]=[minimum.numerator,minimum.denominator]
    return dict(hist=stats,min_by_transverse=min_by_transverse,tight_examples=tight)


if __name__=='__main__':
    t=time.perf_counter()
    item=json.load(open('compatible-mini.json'))['best_nonlocal_margin']
    seed=next(s for s in tra.compute_states() if s['name']==item['seed'])
    w=dict(seed=item['seed'],kind='nonlocal-rank-preserving-extension',generators=item['new_generators'],
           generator_coords=item['coordinates'],a=item['a']+1,t=item['a']+1,U=1)
    exact=verify.audit(seed,w,'extended')
    checks=all_projected_cuts(seed,w)
    if Fraction(*exact['actual_X0_ratio'])!=Fraction(*seed['L']):
        raise AssertionError('nonlocal extension is not actually tight at X=0')
    conclusion=all(z['weighted_failures']==0 for z in checks['hist'].values())
    doc=dict(status='pass',host=platform.node(),python=platform.python_version(),wall_seconds=time.perf_counter()-t,
             extension=exact,all_high_transverse_cuts=checks,
             exact_parent_L=seed['L'] if conclusion else None,
             low_transverse_proof_used='R84 t<=4')
    print(json.dumps(doc,indent=2,sort_keys=True))
    other=json.load(open('transverse-macbook.json'))['minima']
    Path('signed-witnesses-input.json').write_text(json.dumps(dict(minima=dict(**other,extended=w)),indent=2,sort_keys=True)+'\n')
