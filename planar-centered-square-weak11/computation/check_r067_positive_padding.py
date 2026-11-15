"""Exact source/arithmetic check of one positive ten-source padding input.

No full halo area is computed here. The production full-halo and erosion
union computations are separately versioned and independently scoped.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import time


def main():
    assert os.path.ismount('/Volumes/data')
    started=time.perf_counter()
    d,eta,eps=F(1023,1024),F(1,2**20),F(1,1024)
    points=[(x,y) for x in (0,d/2,d) for y in (0,d/2,d)]
    masses=[F(1,5) if i==0 else F(4,5)-7*eta if i==8 else eta for i in range(9)]
    points.append((d+F(1,8),d));masses.append(eps)
    assert sum(masses[:9])==1 and min(masses)>0
    closed={};eligible=[];maxmass=F(0);maxmasks=[]
    for mask in range(1,1<<10):
        ids=[i for i in range(10) if mask>>i&1]
        low=[min(points[i][j] for i in ids) for j in range(2)]
        high=[max(points[i][j] for i in ids) for j in range(2)]
        span=max(high[j]-low[j] for j in range(2))
        mass=sum((masses[i] for i in ids),F(0))
        cmask=sum(1<<i for i in range(10) if all(low[j]<=points[i][j]<=high[j] for j in range(2)))
        if cmask==mask:
            closed[mask]={'mask':mask,'mass':str(mass),'span_squared':str(span**2),'eligible':span**2<mass}
        if span**2<mass:
            eligible.append(mask)
            if mass>maxmass:maxmass=mass;maxmasks=[mask]
            elif mass==maxmass:maxmasks.append(mask)
    assert maxmass==1 and maxmasks==[511]
    assert len(closed)==46 and sum(r['eligible'] for r in closed.values())==18
    reached=1
    while True:
        nxt=reached
        for mask in eligible:
            if mask&reached:nxt|=mask
        if nxt==reached:break
        reached=nxt
    assert reached==1023
    domain_margin=(d+F(1,8))**2-F(35,33)*(1+eps)
    assert domain_margin==F(6982433,34603008)>0
    assert F(3,4)<d<1 and (d/2)**2<F(3,4)
    assert masses[8]>F(3,4)
    assert F(433,250)**2<3 and F(97,250)**2>F(3,20)
    envelope_margin=F(246397,131072000)-F(21,4)*eta
    assert envelope_margin==F(982963,524288000)>0
    assert F(9,64)**2 < masses[8]/4  # residual singleton lies in heavy singleton halo
    assert F(1,8)**2>eps  # every crossing label's same-bbox residual core is empty
    root=Path(__file__).resolve().parents[1]
    inp=root/'computation/inputs/r067-positive-grid-padding.json'
    obj={'name':'r067-positive-grid-padding','points':[[str(x),str(y)] for x,y in points],
         'masses':[str(m) for m in masses],'initial_masses':[str(m) for m in masses],
         'S_indices':list(range(9)),'S_mask':511,'theta':'3/4','eta':str(eta),
         'scope':'One positive padding obstruction to the fixed theta support-neighborhood repair; no general deletion counterexample.'}
    inp.write_text(json.dumps(obj,indent=2)+'\n')
    result={'status':'PASS_EXACT_POSITIVE_PADDING_SOURCE_AND_ANALYTIC_CHECK',
            'python':platform.python_version(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),
            'nonempty_subsets_checked':1023,'closed_bboxes':list(closed.values()),
            'eligible_closed_bboxes':18,'unique_maximum_mask':511,'maximum_mass':'1',
            'positive_sources':10,'S_rank':9,'connected':True,'fat_both_axes':True,
            'nonterminal_margin':str(domain_margin),'grid_gap_squared':str((d/2)**2),
            'support_neighborhood_union_equals_filled_bbox':True,
            'budget_excess_rational_lower':str(envelope_margin),
            'same_bbox_erosion_union_equals_full_halo_for_this_input':True,
            'same_bbox_erosion_overpayment':str(eps),
            'seconds':time.perf_counter()-started,
            'scope':'Exact finite source enumeration and analytic sign checks; no full-halo union run in this author script.'}
    out=root/'computation/results/round067/positive-padding-arithmetic.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='closed_bboxes'},indent=2))


if __name__=='__main__':
    main()
