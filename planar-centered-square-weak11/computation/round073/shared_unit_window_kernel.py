#!/usr/bin/env python3
"""R73 exact finite checks for shared unit-window kernels and non-CDF input.

No area search or Borel n sweep. Uses one n=5 fixture and exact Fraction.
The all-n identities and limits are proved in the accompanying text.
"""
from fractions import Fraction as F
from itertools import combinations,product
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def overlap(a,b,c,d):return max(F(0),min(b,d)-max(a,c))


def kappa(x,w):
    if x<=-1 or x>=w+1:return F(0)
    if x<0:return 1+x
    if x<=w:return F(1)
    return w+1-x


def integrate_piecewise_linear_on_interval(w,a,b):
    """Exact trapezoids at every kink, independently of the closed form."""
    breaks=sorted({a,b,*[x for x in (-F(1),F(0),w,w+1)if a<x<b]})
    return sum(((v-u)*(kappa(u,w)+kappa(v,w))/2 for u,v in zip(breaks,breaks[1:])),F(0))


def closure_catalog(points,masses):
    # Literal coordinate-box membership, independent of previous mask-product
    # enumeration. Every nonempty tight bbox closure has coordinate endpoints.
    xx=sorted({p[0]for p in points});yy=sorted({p[1]for p in points})
    masks=set()
    for i,a in enumerate(xx):
        for b in xx[i:]:
            for j,c in enumerate(yy):
                for d in yy[j:]:
                    masks.add(sum(1<<k for k,p in enumerate(points)if a<=p[0]<=b and c<=p[1]<=d))
    rows=[]
    for mask in sorted(masks-{0}):
        ids=[i for i in range(len(points))if mask>>i&1]
        bbox=[min(points[i][0]for i in ids),max(points[i][0]for i in ids),min(points[i][1]for i in ids),max(points[i][1]for i in ids)]
        mass=sum((masses[i]for i in ids),F(0));span=max(bbox[1]-bbox[0],bbox[3]-bbox[2])
        rows.append({"mask":mask,"ids":ids,"bbox":list(map(str,bbox)),"mass":str(mass),
                     "span_squared":str(span*span),"eligible":mass>span*span,
                     "cap_slack":str(max(F(1),span*span)-mass)})
    return rows


def window(a,b):return(-a,1-a,-b,1-b)


def contains(rect,p):return rect[0]<=p[0]<=rect[1]and rect[2]<=p[1]<=rect[3]


def axis_atoms(edges):
    edges=sorted(set(edges));out=[]
    for i,x in enumerate(edges):
        out.append((x,x,x))
        if i+1<len(edges):out.append((x,edges[i+1],(x+edges[i+1])/2))
    return out


def verify_all_unit_window_deficits(n,e,w,t,h):
    """Whole continuous n=5 unit-center domain, including one-sided limits.

    Counts are constant on point/open interval atoms. Overlap factors are
    affine there, so each bilinear cell extremum lies at a closure corner.
    """
    a,b=-F(1,2),F(3,2);grid=[i*h for i in range(n)]
    edges={a,b}
    edges.update(x+sgn*F(1,2)for x in grid for sgn in (-1,1)if a<=x+sgn*F(1,2)<=b)
    edges.update(x+sgn*F(1,2)for x in (-e,F(1))for sgn in (-1,1)if a<=x+sgn*F(1,2)<=b)
    atoms=axis_atoms(edges);records=[];corner_tests=0;min_g=None;max_g=None
    eps=F(1,n+1)+F(2,n)+2*e
    assert eps<=F(4,n)
    for ix,(x0,x1,xr)in enumerate(atoms):
        k=sum(abs(p-xr)<=F(1,2)for p in grid)
        for x in set((x0,x1)):
            assert abs(F(k,n)-overlap(x-F(1,2),x+F(1,2),0,1))<=F(1,n)
        for iy,(y0,y1,yr)in enumerate(atoms):
            ell=sum(abs(p-yr)<=F(1,2)for p in grid);values=[]
            for x,y in product(set((x0,x1)),set((y0,y1))):
                A=overlap(x-F(1,2),x+F(1,2),-e,1)*overlap(y-F(1,2),y+F(1,2),-e,1)
                g=1-F(k*ell,n*n)-t*(1-A)
                assert 0<=g<=eps
                values.append(g);corner_tests+=1
            low,high=min(values),max(values)
            min_g=low if min_g is None else min(min_g,low)
            max_g=high if max_g is None else max(max_g,high)
            records.append({'axis_atom_indices':[ix,iy],'sigma_coordinate_counts':[k,ell],
                            'deficit_lower_limit':str(low),'deficit_upper_limit':str(high)})
    return {'axis_atoms':[[str(v)for v in a]for a in atoms],
            'all_cell_bounds':records,'cell_count':len(records),'corner_limit_tests':corner_tests,
            'minimum_deficit':str(min_g),'maximum_deficit':str(max_g),
            'proved_uniform_bound':str(eps),'all_centre_boundary_ties_retained':True}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();start=time.perf_counter();raw=args.input.read_bytes();data=json.loads(raw)
    n=data['n'];assert n==5
    e=F(1,n*n);w=1-e;h=F(n+1,n*n);t=F(n,n+1)
    source=[(i*h,j*h)for i in range(n)for j in range(n)]
    q=tuple(map(F,data['outside_point']));outside=F(data['outside_mass'])
    pts=source+[q];masses=[e]*(n*n)+[outside]
    catalog=closure_catalog(pts,masses)
    assert all(F(r['cap_slack'])>=0 for r in catalog)
    eligible=[r for r in catalog if r['eligible']];L=max(F(r['mass'])for r in eligible)
    maxima=[r['mask']for r in eligible if F(r['mass'])==L]
    assert L==1 and maxima==[(1<<(n*n))-1]
    parent=list(range(len(pts)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for r in eligible:
        for i in r['ids'][1:]:parent[find(i)]=find(r['ids'][0])
    components=len({find(i)for i in range(len(pts))});assert components==1
    span=max(max(p[a]for p in pts)-min(p[a]for p in pts)for a in (0,1))
    nonterminal=span*span-F(35,33)*sum(masses,F(0));assert nonterminal>0
    windows=[]
    for key,expected in [('first_window_parameters','first_expected_outside_mass'),('second_window_parameters','second_expected_outside_mass')]:
        a,b=map(F,data[key]);Q=window(a,b);nu=outside if contains(Q,q)else F(0)
        sigma=sum((e for p in source if contains(Q,p)),F(0));assert nu==F(data[expected])
        windows.append({'parameters':[str(a),str(b)],'closed_bounds':list(map(str,Q)),
                        'sigma_mass':str(sigma),'nu_mass':str(nu),'unit_cap_slack':str(1-sigma-nu)})
    assert all(F(a)<=F(b)for a,b in zip(windows[0]['parameters'],windows[1]['parameters']))
    assert F(windows[0]['nu_mass'])>F(windows[1]['nu_mass'])
    C=(-F(1,2),w+F(1,2));J=(-e,F(1))
    # Pointwise exact convolution identity, including every kernel kink and q.
    axis_tests=sorted({-F(3,2),-F(1),-F(1,2),-e,F(0),w,F(1),w+F(1,2),w+1,F(2),q[0],q[1]})
    point_checks=0
    for x,y in product(axis_tests,repeat=2):
        direct=overlap(*C,x-F(1,2),x+F(1,2))*overlap(*C,y-F(1,2),y+F(1,2))
        assert direct==kappa(x,w)*kappa(y,w);point_checks+=1
    internal_kernel_mass=sum((e*kappa(x,w)*kappa(y,w)for x,y in source),F(0));assert internal_kernel_mass==1
    Knu=outside*kappa(q[0],w)*kappa(q[1],w)
    cap_cost=(2-e)**2-1;assert Knu<=cap_cost
    kint_all=integrate_piecewise_linear_on_interval(w,-1,w+1)
    kint_J=integrate_piecewise_linear_on_interval(w,*J)
    assert kint_all==2-e and kint_J==1+e-e*e
    borel_integral=t*(kint_all*kint_all-kint_J*kint_J)
    slack=cap_cost-borel_integral
    formula=((2-e)**2-1+n*((1+e-e*e)**2-1))/F(n+1)
    assert slack==formula and 0<=slack<=F(4,n+1)
    # Four corner-parameter families overcount the central overlap strip.
    corner_cost=3-2*e-e*e
    assert corner_cost-cap_cost==2*e*(1-e)
    all_unit=verify_all_unit_window_deficits(n,e,w,t,h)
    # A very large finite construction is specified by grid indices, not
    # enumerated. These are exact scalar/containment checks for its proof.
    m=100*n*n;delta=F(1,m);lam=1-delta*delta;c=t/(1+2*delta);cprime=lam*c
    assert delta<=F(1,2*n)and delta<=F(1,100)
    forbidden_lo=m-102;forbidden_hi=2*m+1;forbidden_count=forbidden_hi-forbidden_lo+1
    assert forbidden_count==m+104
    assert -1+forbidden_lo*delta==-e-2*delta
    assert -1+(forbidden_hi+1)*delta==1+2*delta
    outside_count=12*m*m-forbidden_count*forbidden_count
    bridge=(-e-F(5,2)*delta,delta/2)
    bridge_distance=max(abs(v)for v in bridge)
    assert bridge_distance<F(1,n)and 2*cprime>1
    coarse_gap=F(399,100)**2-F(455,33)
    assert 4-delta>=F(399,100)and coarse_gap==F(703633,330000)>0
    assert 8*(1+e)*delta+16*delta*delta<=14*delta
    assert t-cprime<=3*delta
    atomic_schematic={'n':n,'mesh_denominator':m,'delta':str(delta),'density_c':str(c),
                      'lambda':str(lam),'density_after_contraction':str(cprime),
                      'B_rectangle':['-1','3','-1','2'],
                      'outside_atom_count_implicit':outside_count,'outside_atoms_enumerated':False,
                      'forbidden_axis_indices':[forbidden_lo,forbidden_hi],
                      'bridge_point':list(map(str,bridge)),
                      'bridge_squared_distance_gap_to_internal_mass':str(e-bridge_distance*bridge_distance),
                      'neighbor_pair_eligibility_factor_margin':str(2*cprime-1),
                      'nonterminal_uniform_lower_gap':str(coarse_gap),
                      'all_unit_window_mass_comparison_error_bound':str(21*delta),
                      'caps_source':'separate R73 cap-preserving atomization theorem; not inferred from these scalar checks'}
    out={'schema':'r073-shared-unit-window-kernel-result-v1','all_passed':True,
         'n':n,'source_count':len(pts),'closure_count':len(catalog),'eligible_closure_count':len(eligible),
         'maximum_mass':str(L),'maximum_masks':maxima,'eligible_component_count':components,
         'nonterminal_gap':str(nonterminal),'windows':windows,
         'cdf_monotonicity_violated_for_actual_admissible_input':True,
         'pointwise_kernel_checks':point_checks,'internal_kernel_mass':str(internal_kernel_mass),
         'uniform_C_kernel_budget':str(cap_cost),'one_atom_weighted_kernel_integral':str(Knu),
         'Borel_family_kernel_integral':str(borel_integral),'Borel_integrated_unit_cap_slack':str(slack),
         'all_n_formula':'[(3-4e+e^2)+n((1+e-e^2)^2-1)]/(n+1), e=1/n^2',
         'all_n_upper_bound':'4/(n+1); all-n proof in accompanying note, not a finite scan',
         'symmetric_four_corner_budget':str(corner_cost),'corner_overcount_cost':str(corner_cost-cap_cost),
         'complete_Borel_unit_window_center_certificate':all_unit,
         'finite_atomic_schematic_checks':atomic_schematic,
         'complete_closure_catalog':catalog,
         'actual_atoms':{'threshold':'1','positions':[[str(v)for v in p]for p in pts],
                         'masses':list(map(str,masses)),'S_indices':list(range(n*n))},
         'scope':'exact window/kernel bookkeeping and valid finite non-CDF counterexample; no complete halo area computed, no 3.5/C2 theorem; Borel slack obstruction not promoted to finite-atomic area claim',
         'receipt':{'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),
                    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start}}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['all_passed','source_count','closure_count','eligible_closure_count','maximum_masks','nonterminal_gap','uniform_C_kernel_budget','Borel_integrated_unit_cap_slack']}))


if __name__=='__main__':main()
