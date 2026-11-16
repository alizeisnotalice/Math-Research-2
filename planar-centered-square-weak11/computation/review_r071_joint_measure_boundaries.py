#!/usr/bin/env python3
"""Literal Fraction controls for R71; no producer LP/geometry imports."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time


def dist(p, q):
    return max(abs(p[0]-q[0]), abs(p[1]-q[1]))


def bbox(points):
    return min(x for x,y in points), max(x for x,y in points), min(y for x,y in points), max(y for x,y in points)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True); args=ap.parse_args(); start=time.perf_counter()
    points=[(F(0),F(0)),(F(1,2),F(0)),(F(0),F(3,2)),(F(1,2),F(3,2))]
    receivers=[(F(0),F(0)),(F(1,4),F(0)),(F(1,4),F(3,4)),(F(2),F(0)),(F(1,2),F(3,2)),(F(-1,4),F(3,4))]
    weights=[F(0),F(1,4),F(1,2),F(1)]
    masks=[tuple(i for i in range(len(points)) if mask>>i&1) for mask in range(1,1<<len(points))]
    closures={}
    for subset in masks:
        a,b,c,d=bbox([points[i] for i in subset]); mask=tuple(i for i,(x,y) in enumerate(points) if a<=x<=b and c<=y<=d)
        closures[mask]=max(b-a,d-c)
    stats=dict(mass_assignments=0,literal_subset_eligibility_checks=0,bbox_cap_rows=0,receiver_cases=0,literal_halo_subset_checks=0,complete_tied_prefix_rows=0,positive_zero_distance_cases=0)
    for masses in product(weights,repeat=len(points)):
        stats['mass_assignments']+=1
        largest=F(0)
        for subset in masks:
            positive=[i for i in subset if masses[i]>0]
            if not positive: continue
            a,b,c,d=bbox([points[i] for i in positive]); mass=sum((masses[i] for i in positive),F(0)); diameter=max(b-a,d-c)
            stats['literal_subset_eligibility_checks']+=1
            if diameter**2<mass: largest=max(largest,mass)
        cap=True
        for subset,diameter in closures.items():
            stats['bbox_cap_rows']+=1
            cap &= sum((masses[i] for i in subset),F(0))<=max(F(1),diameter**2)
        assert cap==(largest<=1),(masses,cap,largest)
        for x in receivers:
            stats['receiver_cases']+=1
            direct=False
            for subset in masks:
                positive=[i for i in subset if masses[i]>0]
                if not positive: continue
                mass=sum((masses[i] for i in positive),F(0)); R=max(dist(points[i],x) for i in positive)
                stats['literal_halo_subset_checks']+=1
                direct |= 4*R**2<mass
            prefix=False
            for R in sorted(set(dist(q,x) for q in points)):
                stats['complete_tied_prefix_rows']+=1
                mass=sum((masses[i] for i,q in enumerate(points) if dist(q,x)<=R),F(0)); prefix |= mass>4*R**2
                if R==0 and mass>0: stats['positive_zero_distance_cases']+=1
            assert direct==prefix,(masses,x,direct,prefix)
    # Boundary capture: closed radius captures mass but open radius does not.
    r=F(1,4); mass=F(1,2); R=F(1,4); perturbed=F(3,10)
    assert R==r and not R<r and mass>4*r*r
    assert R<perturbed and mass>4*perturbed**2 and perturbed<F(1,2)
    # An artificial algebraic shared-cap example, not a geometric realization.
    y=F(1,2); lambdas=[F(1,2),F(1,2)]; demands=[F(3,5),F(3,5)]
    assert sum(lambdas)==1 and y>=lambdas[0] and y>=lambdas[1]
    dual=y-sum((l*d for l,d in zip(lambdas,demands)),F(0)); assert dual==F(-1,10)
    assert F(7,10)<=1 and F(7,10)-F(3,5)==F(1,10) # each singleton has eta=.1
    result=dict(status='PASS_LITERAL_R71_BOUNDARY_EQUIVALENCE_CONTROLS',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version,platform=platform.platform(),stats=stats,positions=[[str(c) for c in q] for q in points],receivers=[[str(c) for c in q] for q in receivers],mass_values=[str(w) for w in weights],bbox_closures=len(closures),boundary_example=dict(radius=str(r),source_distance=str(R),source_mass=str(mass),larger_true_witness=str(perturbed)),algebraic_toy_dual_eta_upper=str(dual),elapsed_seconds=time.perf_counter()-start,scope='Only exact finite equivalence/boundary controls and algebraic dual signs; no geometric shared-measure obstruction, LP production replay, halo area or new C2 bound.')
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(dict(status=result['status'],stats=stats,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__': main()
