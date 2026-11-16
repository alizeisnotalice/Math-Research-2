#!/usr/bin/env python3
"""Exact actual-atom whole-cell certificate check, independent of MIP rows."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
from review_r072_actual_halo_y_slabs import direct_closures,y_slabs


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def delta(bb,C):
    a,b,c,d=bb;u,v,s,t=C
    return max((b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,F(0))
def Rmax(q,C):
    u,v,s,t=C
    return max(abs(q[0]-u),abs(q[0]-v),abs(q[1]-s),abs(q[1]-t))


def main():
    ap=argparse.ArgumentParser()
    for n in ('model','result','actual_input','contracted_input','output'):ap.add_argument('--'+n.replace('_','-'),type=Path,required=True)
    a=ap.parse_args();start=time.perf_counter();m=json.loads(a.model.read_text());saved=json.loads(a.result.read_text());cert=saved['certificate'];raw=json.loads(a.actual_input.read_text());contracted=json.loads(a.contracted_input.read_text());den=m['units_denominator']
    S=[tuple(F(c,den) for c in q) for q in m['S_positions_units']];O=[tuple(F(c,den) for c in q) for q in m['outside_positions_units']];weights=list(map(F,cert['outside_weights']));assert len(S)==25 and len(O)==len(weights)==145 and min(weights)>=0
    assert len(set(S+O))==len(S)+len(O)
    pos=[q for q,w in zip(O,weights) if w>0];pm=[w for w in weights if w>0]
    expected=S+pos;em=[F(1,25)]*25+pm
    assert [tuple(map(F,q)) for q in raw['positions']]==expected and list(map(F,raw['masses']))==em
    assert [tuple(map(F,q)) for q in contracted['positions']]==expected
    cm=list(map(F,contracted['masses']));assert cm[:25]==em[:25]
    ratios={b/a for a,b in zip(em[25:],cm[25:])};assert ratios=={F(99999999,100000000)};lam=ratios.pop()
    residual,scans=direct_closures(pos,pm);full,fs=direct_closures(expected,em)
    L=max(r['mass'] for r in full if r['eligible']);assert L==1
    Cfull=[tuple(F(c,den) for c in box) for box in m['cells_units']]
    (partition_area,_,_),_=y_slabs(Cfull,[]);sum_partition=sum((b-a)*(d-c) for a,b,c,d in Cfull);assert partition_area==sum_partition
    # Pairwise positive-area intersections are forbidden, shared edges allowed.
    checks=0
    for i,C in enumerate(Cfull):
        for D in Cfull[:i]:
            assert min(C[1],D[1])<=max(C[0],D[0]) or min(C[3],D[3])<=max(C[2],D[2]);checks+=1
    rows=[];selected=[];total_residual_checks=0;minmargin=None;contractmin=None
    for claim in cert['selected_cells']:
        C=tuple(map(F,claim['bounds']));assert C==Cfull[claim['cell_id']]
        assert len(set(claim['label_S_indices']))==len(claim['label_S_indices']) and len(set(claim['label_outside_indices']))==len(claim['label_outside_indices'])
        label=[S[i] for i in claim['label_S_indices']]+[O[i] for i in claim['label_outside_indices']]
        assert label
        mass=F(len(claim['label_S_indices']),25)+sum((weights[i] for i in claim['label_outside_indices']),F(0))
        r=max(Rmax(q,C) for q in label);assert r==F(claim['closed_full_label_radius'])<F(1,2)
        margin=mass-4*r*r;assert mass==F(claim['captured_mass']) and margin==F(claim['strict_full_margin'])>=F(m['strict_full_margin'])>0
        c_mass=F(len(claim['label_S_indices']),25)+lam*sum((weights[i] for i in claim['label_outside_indices']),F(0));c_margin=c_mass-4*r*r;assert c_margin>0
        worst=None
        for T in residual:
            slack=4*delta(T['bbox'],C)**2-T['mass'];assert slack>=0
            worst=slack if worst is None else min(worst,slack);total_residual_checks+=1
        # All residual constraints improve after outside contraction.
        minmargin=margin if minmargin is None else min(minmargin,margin);contractmin=c_margin if contractmin is None else min(contractmin,c_margin)
        selected.append(C);rows.append(dict(cell_id=claim['cell_id'],mass=str(mass),whole_cell_max_label_radius=str(r),strict_full_margin=str(margin),contracted_strict_full_margin=str(c_margin),minimum_residual_closure_slack=str(worst)))
    (unionarea,_,_),stats=y_slabs(selected,[]);sumarea=sum((b-a)*(d-c) for a,b,c,d in selected)
    assert unionarea==sumarea==F(cert['certified_area'])==F(345909,250000)
    assert len(selected)==cert['selected_count']==89 and minmargin==F(cert['smallest_strict_full_margin'])
    result=dict(status='PASS_ACTUAL_POSITIVE_ATOM_WHOLE_CELL_DELETION_CERTIFICATE',source_sha256=sha(Path(__file__)),independent_closure_and_area_helper_sha256=sha(Path(__file__).with_name('review_r072_actual_halo_y_slabs.py')),input_sha256={n:sha(getattr(a,n)) for n in ('model','result','actual_input','contracted_input')},python=sys.version,actual_positive_sources=len(expected),actual_positive_outside=len(pos),maximum_eligible_mass=str(L),actual_full_closures=len(full),actual_residual_closures=len(residual),actual_residual_coordinate_rectangles=scans,all_cells=len(Cfull),partition_pair_checks=checks,selected_closed_cells=len(selected),literal_selected_cell_residual_closure_checks=total_residual_checks,minimum_strict_full_margin=str(minmargin),contracted_minimum_strict_full_margin=str(contractmin),contracted_weight_factor=str(lam),certified_actual_Delta_area_lower=str(unionarea),area_slab_stats=stats,cells=rows,elapsed_seconds=time.perf_counter()-start,scope='Actual positive atoms and every actual residual bbox closure at every selected cell; strict full-label whole-cell distances including zero candidate label positions; independent union area. Same cells remain valid under contracted weights. No MIP optimum/dual/menu-completeness or general upper bound claimed.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],cells=len(selected),actual_residual_closures=len(residual),all_residual_cell_checks=total_residual_checks,area=str(unionarea),seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
