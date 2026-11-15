#!/usr/bin/env python3
"""R64 exact independent constraints + 2D compressed-grid deletion audit.

No producer imports, finite differences, optimizer, or MATLAB. Geometry is
rebuilt for seed and best only; all eight recorded line masses are checked.
"""
from fractions import Fraction as F
import hashlib,json,os
from itertools import combinations_with_replacement
from pathlib import Path
import subprocess,sys,time

ROOT=Path('/Volumes/data/CodexProjects/T2')
RUN=ROOT/'.local/compute/r064-lp-pilot'
AUTHOR=ROOT/'.local/workspaces/pcw-astra-covering/planar-centered-square-weak11'
PRODUCER=ROOT/'.local/workspaces/pcw-position-slab-impl/planar-centered-square-weak11/computation/round064/mass_polytope_variation.py'
HERE=Path(__file__).resolve().parent
CODE='4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246'
INPUT='46acda723e8590054912a18f71c67efc95b5a0bcfd449a9c60ce615413ac41dd'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)


def boxes(points):
    result={}
    xs=sorted({x for x,y in points});ys=sorted({y for x,y in points})
    for a,b in combinations_with_replacement(xs,2):
        for c,d in combinations_with_replacement(ys,2):
            ids=tuple(i for i,(x,y)in enumerate(points)if a<=x<=b and c<=y<=d)
            if not ids:continue
            tight=(min(points[i][0]for i in ids),max(points[i][0]for i in ids),min(points[i][1]for i in ids),max(points[i][1]for i in ids))
            if tight!=(a,b,c,d):continue
            mask=sum(1<<i for i in ids);assert mask not in result
            result[mask]=dict(ids=ids,bbox=tight,d2=max(b-a,d-c)**2)
    return result


def span2(points,ids):
    return max(max(points[i][a]for i in ids)-min(points[i][a]for i in ids)for a in (0,1))**2


def feasible(points,mass,spec,bb):
    S=spec['S_mask'];ss=spec['S_indices'];C=spec['mass_constraints'];eta=F(spec['parameters']['unique_margin_eta'])
    assert S==sum(1<<i for i in ss) and sum(mass[i]for i in ss)==1 and span2(points,ss)<1
    floor=F(C['positive_floor_for_initial_chamber']);assert min(mass)>=floor>0
    eligible=[];sums={}
    for mask,b in bb.items():
        m=sum(mass[i]for i in b['ids']);sums[mask]=m
        assert m<=max(F(1),b['d2'])
        if mask!=S and b['d2']<1:assert m<=1-eta
        if m>b['d2']:eligible.append(mask)
    maximum=max(sums[s]for s in eligible)
    winners=sorted(s for s in eligible if sums[s]==maximum)
    assert maximum==1 and winners==[S]
    # Actual hypergraph connectivity, independently recomputed by union-find.
    parent=list(range(len(points)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for mask in eligible:
        ids=bb[mask]['ids'];a=find(ids[0])
        for i in ids[1:]:parent[find(i)]=a
    assert len({find(i)for i in range(len(points))})==1
    total=sum(mass);d2=span2(points,range(len(points)))
    cap=F(C['nonterminal_total_mass_cap']);assert cap==F(33,35)*(d2-F(1,16)) and total<=cap
    margin=d2-F(35,33)*total;assert margin>=F(1,16)
    conn=F(C['connectivity_margin']);covered=set()
    for edge in C['prescribed_connection_lower_constraints']:
        ids=[i for i in range(len(points))if edge['mask']&(1<<i)]
        ed2=span2(points,ids)
        assert F(edge['bbox_span_squared'])==ed2
        assert F(edge['minimum_mass'])==ed2+conn
        assert sum(mass[i]for i in ids)>=ed2+conn
        assert edge['outside_index']in ids and edge['mask']&S
        covered.add(edge['outside_index'])
    assert covered==set(range(len(points)))-set(ss)
    a,b,c,d=bb[S]['bbox']
    for i,p in enumerate(points):
        assert max(F(0),a-p[0],p[0]-b,c-p[1],p[1]-d)<=1
        if not S&(1<<i):assert max(b-p[0],p[0]-a,d-p[1],p[1]-c)>1
    return dict(positive_sources=len(mass),closures=len(bb),eligible_closures=len(eligible),maximum_mass='1',winners=winners,
                actual_connected=True,actual_nonterminal=True,total_mass=str(total),actual_span2=str(d2),nonterminal_margin=str(margin),
                all_bbox_caps=True,unique_gap=True,positive_floor=True,prescribed_edges_verified=len(covered))


def root_floor(num,den):
    lo,hi=0,1<<((num.bit_length()+1)//2)
    while lo+1<hi:
        mid=(lo+hi)//2
        if den*mid*mid<=num:lo=mid
        else:hi=mid
    assert den*lo*lo<=num<den*(lo+1)*(lo+1)
    return lo


def rectangles(points,mass,bb,bits):
    D=1<<bits;out=[[],[]]
    for b in bb.values():
        m=sum(mass[i]for i in b['ids'])
        if m<=b['d2']:continue
        k=root_floor(m.numerator*D*D,4*m.denominator)
        a,bx,c,d=[v*D for v in b['bbox']];assert all(q.denominator==1 for q in (a,bx,c,d))
        a,bx,c,d=map(int,(a,bx,c,d))
        for which,r in ((0,max(0,k-1)),(1,k+1)):
            rec=(bx-r,a+r,d-r,c+r)
            if rec[0]<rec[1]and rec[2]<rec[3]:out[which].append(rec)
    return out


def compressed_area(layers):
    """Four independent integer 2D difference arrays; no interval merging."""
    xs=sorted({r[j]for layer in layers for r in layer for j in (0,1)})
    ys=sorted({r[j]for layer in layers for r in layer for j in (2,3)})
    nx,ny=len(xs),len(ys);xi={x:i for i,x in enumerate(xs)};yi={y:i for i,y in enumerate(ys)}
    arrays=[[0]*(nx*ny)for _ in layers]
    for ar,layer in zip(arrays,layers):
        for xl,xr,yl,yr in layer:
            l,r,a,b=xi[xl],xi[xr],yi[yl],yi[yr]
            for x,y,v in ((l,a,1),(r,a,-1),(l,b,-1),(r,b,1)):ar[x*ny+y]+=v
    areas=[0]*6;cells=0
    for ix in range(nx):
        for iy in range(ny):
            j=ix*ny+iy
            for ar in arrays:
                if ix:ar[j]+=ar[j-ny]
                if iy:ar[j]+=ar[j-1]
                if ix and iy:ar[j]-=ar[j-ny-1]
                assert ar[j]>=0
            if ix+1==nx or iy+1==ny:continue
            covered=[ar[j]>0 for ar in arrays]
            assert not covered[0]or covered[1]
            assert not covered[2]or covered[3]
            area=(xs[ix+1]-xs[ix])*(ys[iy+1]-ys[iy]);cells+=1
            for k in range(4):
                if covered[k]:areas[k]+=area
            if covered[0]and not covered[3]:areas[4]+=area
            if covered[1]and not covered[2]:areas[5]+=area
    return areas,cells


def objective(points,mass,bb,S,bits):
    full=rectangles(points,mass,bb,bits)
    keep=[i for i in range(len(points))if not S&(1<<i)]
    pp=[points[i]for i in keep];mm=[mass[i]for i in keep];rb=boxes(pp)
    residual=rectangles(pp,mm,rb,bits)
    ar,cells=compressed_area(full+residual);D2=1<<(2*bits)
    result={k:str(F(v,D2))for k,v in zip(('full_inner_area','full_outer_area','remaining_inner_area','remaining_outer_area','lower','upper'),ar)}
    result.update(full_active_boxes=len(full[1]),remaining_active_boxes=len(residual[1]),bits=bits,compressed_grid_cells=cells)
    # Exact similarity check: scale positions16/masses256 and reduce bits4.
    scaled=[(16*x,16*y)for x,y in points];scaledmass=[256*m for m in mass]
    sb=boxes(scaled);sf=rectangles(scaled,scaledmass,sb,bits-4)
    sp=[scaled[i]for i in keep];sm=[scaledmass[i]for i in keep]
    sr=rectangles(sp,sm,boxes(sp),bits-4)
    assert [sorted(x)for x in full+residual]==[sorted(x)for x in sf+sr]
    result['scale16_mass256_exact_integer_rectangles_match']=True
    return result


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    assert sha(PRODUCER)==sha(RUN/'mass_polytope_variation.py')==CODE
    assert sha(RUN/'input.json')==sha(AUTHOR/'computation/inputs/r064-fixed-position-positive21.json')==INPUT
    result=json.loads((RUN/'result.json').read_text());receipt=json.loads((RUN/'receipt.json').read_text())
    assert result['source_sha256']==receipt['source_sha256']==CODE
    assert result['input_sha256']==receipt['input_sha256']==INPUT
    assert receipt['exit_status']==0 and receipt['owned_pid_absent']is True
    assert not subprocess.run(['ps','-p',str(receipt['pid']),'-o','pid='],capture_output=True,text=True).stdout.strip()
    spec=json.loads((RUN/'input.json').read_text());points=[tuple(map(F,p))for p in spec['points']]
    seed=list(map(F,spec['initial_masses']));bb=boxes(points);S=spec['S_mask']
    seedcheck=feasible(points,seed,spec,bb);assert len(bb)==274 and seedcheck['eligible_closures']==54
    linechecks=[];direction=list(map(F,result['rational_direction']))
    for record in result['candidates']:
        mass=list(map(F,record['masses']));step=F(record['step'])
        assert mass==[m+step*d for m,d in zip(seed,direction)]
        checked=feasible(points,mass,spec,bb)
        assert checked['eligible_closures']==record['exact_check']['strict_eligible_closures']
        assert record['exact_check']['feasible']and not record['exact_check']['errors']
        assert record['certified_strict_improvement']is True
        assert F(record['certified_gain_lower'])==F(record['objective']['lower'])-F(result['seed_objective']['upper'])>0
        linechecks.append(dict(step=str(step),exact_constraint_check=checked,reported_gain_arithmetic_verified=True,
                               geometry_independently_recomputed=record['step']==result['best_candidate']['step']))
    best=result['best_candidate'];bestmass=list(map(F,best['masses']));assert best['step']=='1'
    bits=result['parameters']['bits'];indseed=objective(points,seed,bb,S,bits);indbest=objective(points,bestmass,bb,S,bits)
    for independent,saved in ((indseed,result['seed_objective']),(indbest,best['objective'])):
        for key,value in saved.items():
            if key!='scope':assert independent[key]==value,(key,independent[key],value)
    gain=F(indbest['lower'])-F(indseed['upper']);assert gain>0 and str(gain)==best['certified_gain_lower']
    proof=AUTHOR/'proofs/r064-source-state-mass-polytope-variation.md'
    out=dict(status='PASS_R64_PROOF_GEOMETRY_AND_BEST_FINITE_STEP',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,
             source_pins=[pin(p)for p in (PRODUCER,RUN/'input.json',RUN/'result.json',RUN/'receipt.json',proof)],
             arithmetic='Exact Fraction/int; coordinate-tight bbox enumeration, integer bisection roots, 4-layer compressed2D difference-array area. No producer import, LP/gradient execution, MATLAB or parameter search.',
             seed_constraints=seedcheck,all_eight_line_mass_constraints=linechecks,seed_geometry=indseed,best_geometry=indbest,
             independently_certified_gain_lower=str(gain),root_receipt_exit_status=0,root_pid_absent_now=True,
             scope='All8 mass/design/support/connectivity/nonterminal checks independent; complete geometric areas independently recomputed only seed+best. Other7 line geometric certificates not replayed. Crossing eligibility walls accepted as finite-step improvement, no chamber-gradient/local/global optimality/general C2 claim.')
    target=HERE/'results/round064-mass-variation-independent.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('all_eight_line_mass_constraints','source_pins')},indent=2))


if __name__=='__main__':main()
