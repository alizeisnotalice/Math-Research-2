#!/usr/bin/env python3
"""Independent R53 whole simplex BSP and inclusion-exclusion area audit."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, platform, sys, time

def area_inclusion_exclusion(rectangles):
    value=0
    for mask in range(1,1<<len(rectangles)):
        active=[rectangles[i] for i in range(len(rectangles)) if mask>>i&1]
        dx=max(0,min(r[1] for r in active)-max(r[0] for r in active))
        dy=max(0,min(r[3] for r in active)-max(r[2] for r in active))
        value += (1 if mask.bit_count()%2 else -1)*dx*dy
    return value

def run(cert):
    records=cert['nodes'];D=cert['coordinate_denominator'];assert D==2**cert['bits'] and D%4==0
    bands=((Q(0),Q(1,4)),(Q(1,4),Q(3,4)),(Q(3,4),Q(1)))
    seen=set();leaves=[];empty=0;depthmax=0
    def visit(index,oldlo,oldhi,depth):
        nonlocal empty,depthmax
        assert isinstance(index,int) and 0<=index<len(records) and index not in seen
        seen.add(index);record=records[index]
        # Every removed interval is impossible under sum(w)=1. One simultaneous pass.
        lo=[max(oldlo[j],1-sum(oldhi[k] for k in range(9) if k!=j)) for j in range(9)]
        hi=[min(oldhi[j],1-sum(oldlo[k] for k in range(9) if k!=j)) for j in range(9)]
        feasible=sum(lo)<=1<=sum(hi) and all(l<=h for l,h in zip(lo,hi))
        if record['kind']=='empty':
            assert not feasible;empty+=1;return
        assert feasible and list(map(Q,record['lo']))==lo and list(map(Q,record['hi']))==hi
        assert all(0<=l<=h<=1 for l,h in zip(lo,hi))
        if record['kind']=='split':
            j=record['coordinate'];c=Q(record['cut']);assert 0<=j<9 and lo[j]<c<hi[j]
            # Closed children overlap at the cut, so even boundary faces are covered.
            hi_left=hi.copy();hi_left[j]=c
            lo_right=lo.copy();lo_right[j]=c
            visit(record['left'],lo,hi_left,depth+1)
            visit(record['right'],lo_right,hi,depth+1)
            return
        assert record['kind']=='leaf' and record['depth']==depth
        order=record['order'];assert sorted(order)==list(range(9))
        assert len(record['outer_rectangles_integer'])==9
        rectangles=[]
        for k,j in enumerate(order):
            prefix=set(order[:k+1]);suffix=set(range(9))-prefix
            upper=min(sum(hi[t] for t in prefix),1-sum(lo[t] for t in suffix))
            assert 0<=upper<=1 and upper==Q(record['prefix_upper'][k])
            x0,x1=bands[j//3];y0,y1=bands[j%3]
            rr=tuple(record['outer_rectangles_integer'][k]);assert all(isinstance(x,int) for x in rr)
            radius=int(x0*D)-rr[0]
            assert radius>=0 and rr==(int(x0*D)-radius,int(x1*D)+radius,int(y0*D)-radius,int(y1*D)+radius)
            # Direct integer square inequality, without generator isqrt or ceil routine.
            assert 4*radius*radius*upper.denominator>=D*D*upper.numerator
            assert radius==0 or 4*(radius-1)**2*upper.denominator<D*D*upper.numerator
            rectangles.append(rr)
        numerator=area_inclusion_exclusion(rectangles)
        assert numerator==record['area_upper_numerator'] and record['area_denominator']==D*D
        assert 2*numerator<7*D*D
        leaves.append((index,Q(numerator,D*D)));depthmax=max(depthmax,depth)
    visit(0,[Q(0)]*9,[Q(1)]*9,0)
    assert seen==set(range(len(records)))
    assert len(records)==cert['node_count'] and len(leaves)==cert['leaf_count'] and depthmax==cert['max_depth']
    bound=max(v for _,v in leaves)
    assert bound==Q(cert['uniform_area_upper_exact']) and Q(7,2)-bound==Q(cert['strict_margin_below_7_2'])
    assert dict(leaves)[cert['maximum_leaf']]==bound
    # Exact scalar branch checks for the separately reviewed contraction error bound.
    scalar=0
    for rootmass in [Q(1,100),Q(1,4),Q(1),Q(7,2)]:
      for delta in [Q(0),Q(1,100),Q(1,10),Q(1)]:
       for sx in [Q(0),rootmass/2,rootmass,rootmass+delta,rootmass+2*delta,rootmass+3*delta]:
        for sy in [Q(0),rootmass/2,rootmass,rootmass+delta,rootmass+2*delta,rootmass+3*delta]:
         e=max(0,rootmass+2*delta-sx)*max(0,rootmass+2*delta-sy)-max(0,rootmass-sx)*max(0,rootmass-sy)
         assert 0<=e<=4*delta*rootmass+4*delta*delta;scalar+=1
    return dict(status='PASS',nodes=len(records),leaves=len(leaves),empty_nodes=empty,maximum_depth=depthmax,union_method='Exact inclusion-exclusion over all 511 nonempty subsets of the nine outer rectangles, independent of author sweeps and grids.',rectangle_subsets_checked=511*len(leaves),uniform_area_upper_exact=str(bound),strict_margin_below_7_2=str(Q(7,2)-bound),root_domain='[0,1]^9 intersect sum(w)=1; closed branches include all faces and zero masses',outward_radius_checks=9*len(leaves),contraction_scalar_checks=scalar)

def main():
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();t=time.perf_counter();r=run(json.loads(a.certificate.read_text()))
    r.update(python=sys.version,platform=platform.platform(),elapsed_seconds=time.perf_counter()-t,source_sha256={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.certificate,Path(__file__).resolve()]},scope='Whole finite mass-simplex certificate. General arbitrary-bbox C2 remains open; continuous theorem mapping and contraction reviewed analytically.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
