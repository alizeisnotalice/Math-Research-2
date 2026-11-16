#!/usr/bin/env python3
"""Independent exact source/cap/prefix/dual audit. No producer imports or LP."""
import argparse
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from pathlib import Path
import hashlib
import json
import sys
import time


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def inside(p,b): return b[0]<=p[0]<=b[1] and b[2]<=p[1]<=b[3]
def dist(p,x): return max(abs(p[0]-x[0]),abs(p[1]-x[1]))
def box(x,r): return x[0]-r,x[0]+r,x[1]-r,x[1]+r
def fs(p): return [str(v) for v in p]


def all_caps(S,O):
    """All coordinate rectangles, min cap for each outside capture mask."""
    axes=[sorted({p[k] for p in S+O}) for k in (0,1)]
    result={}; count=0
    for (a,b),(c,d) in product(combinations_with_replacement(axes[0],2),combinations_with_replacement(axes[1],2)):
        bb=(a,b,c,d); hit=tuple(i for i,q in enumerate(O) if inside(q,bb))
        if not hit: continue
        count+=1
        old=sum(inside(p,bb) for p in S)*F(1,25)
        cap=max(F(1),max(b-a,d-c)**2)-old
        if hit not in result or cap<result[hit]: result[hit]=cap
    return result,count


def load_model(path,S,O,centres,radius,original):
    obj=json.loads(path.read_text()); rows=[]; metas=[]; seen_caps={}; seen_residual={x:set() for x in centres}; witnesses=[]
    for rr,vv,label,index in zip(obj['rows'],obj['values'],obj['labels'],range(len(obj['rows']))):
        kind=label[0]
        if original:
            if kind=='cap':
                _,il,ih,jl,jh,ns=label; bb=tuple(F(v*24,100) for v in (il,ih,jl,jh)); assert ns==sum(inside(p,bb) for p in S)
            else:
                _,x,y,r,*extra=label; x=(F(x,100),F(y,100)); r=F(r,100); bb=box(x,r)
                if extra: assert extra==[sum(inside(p,bb) for p in S)]
            rhs=F(obj['rhs_numerators'][index],obj['rhs_denominator'])
        else:
            bb=tuple(map(F,label[1:])); rhs=F(obj['rhs'][index])
        hit=tuple(i for i,q in enumerate(O) if inside(q,bb))
        if kind in ('cap','collar'):
            cap=max(F(1),max(bb[1]-bb[0],bb[3]-bb[2])**2)-sum(inside(p,bb) for p in S)*F(1,25)
            assert rhs==cap and rr==list(hit) and vv==[1]*len(hit)
            if kind=='cap': seen_caps[hit]=min(seen_caps.get(hit,cap),cap)
            else: assert cap==0
        else:
            x=((bb[0]+bb[1])/2,(bb[2]+bb[3])/2); r=(bb[1]-bb[0])/2
            assert bb[3]-bb[2]==2*r and x in centres
            if kind=='residual':
                assert rhs==4*r*r and rr==list(hit) and vv==[1]*len(hit)
                seen_residual[x].add(r)
            elif kind=='witness':
                assert r==radius and rhs==sum(inside(p,bb) for p in S)*F(1,25)-4*r*r
                assert rr==list(hit)+[len(O)] and vv==[-1]*len(hit)+[1]
                witnesses.append(x)
            else: raise AssertionError(kind)
        rows.append((rr,vv,rhs)); metas.append((kind,bb))
    assert sorted(witnesses)==sorted(centres)
    for x in centres: assert seen_residual[x]=={dist(q,x) for q in O}
    expected,count=all_caps(S,O); assert seen_caps==expected,(len(seen_caps),len(expected))
    return rows,metas,dict(actual_export_rows=len(rows),all_coordinate_rectangles_with_outside_mass=count,distinct_minimal_outside_cap_masks=len(expected),all_complete_residual_prefixes=sum(len(t) for t in seen_residual.values()),strict_witness_rows=len(witnesses))


def verify_lp(rows,cert,n):
    z=list(map(F,cert['primal']['weights'])); eta=F(cert['primal']['eta']); assert len(z)==n and min(z)>=0
    allz=z+[eta]; slacks=[]
    for rr,vv,b in rows:
        slack=b-sum((allz[j]*v for j,v in zip(rr,vv)),F(0)); assert slack>=0; slacks.append(slack)
    weights=[F(0)]*len(rows)
    for i,a in cert['dual']['nonzero_rows']: weights[i]+=F(a)
    assert min(weights)>=0
    cols=[F(0)]*(n+1)
    bound=F(0)
    for (rr,vv,b),a in zip(rows,weights):
        if not a: continue
        bound+=a*b
        for j,v in zip(rr,vv): cols[j]+=a*v
    assert cols[n]==1 and min(cols[:n])>=0 and bound==eta
    assert bound==F(cert['dual'].get('upper_eta', cert['dual'].get('upper_eta')))
    return dict(eta=str(eta),dual_bound=str(bound),positive_outside_sources=sum(a>0 for a in z),all_rows_primal_feasible=True,all_dual_columns_nonnegative=True,free_eta_dual_column='1',minimum_row_slack=str(min(slacks)))


def actual_support(S,O,z,centres,radius,strict):
    pos=[q for q,a in zip(O,z) if a>0]; weights=[a for a in z if a>0]
    caps,count=all_caps(S,pos)
    for hit,cap in caps.items(): assert sum((weights[j] for j in hit),F(0))<=cap
    answer=[]
    for x in centres:
        residual=[]
        for R in sorted({dist(q,x) for q in pos}):
            M=sum((a for q,a in zip(pos,weights) if dist(q,x)<=R),F(0)); assert M<=4*R*R
            residual.append((R,4*R*R-M))
        old_prefixes=[(R,sum(dist(p,x)<=R for p in S)*F(1,25)) for R in sorted({dist(p,x) for p in S})]
        assert all(M<=4*R*R for R,M in old_prefixes) # these receivers are outside the old halo too
        captured_old=[p for p in S if dist(p,x)<=radius]
        captured_ext=[(q,a) for q,a in zip(pos,weights) if dist(q,x)<=radius]
        mass=len(captured_old)*F(1,25)+sum((a for _,a in captured_ext),F(0)); margin=mass-4*radius**2
        row=dict(receiver=fs(x),residual_distance_events=len(residual),residual_nonhalo=True,old_nonhalo=True,closed_witness_mass=str(mass),closed_witness_margin=str(margin))
        if strict:
            assert margin>0 and captured_old and captured_ext
            open_radius=radius+min((F(1,2)-radius)/2,margin/16)
            assert radius<open_radius<F(1,2) and mass>4*open_radius**2
            row['genuine_larger_open_witness_radius']=str(open_radius)
            row['strict_positive_margin']=True
        answer.append(row)
    return dict(actual_positive_outside=len(pos),actual_positive_total=len(S)+len(pos),actual_coordinate_rectangles_with_outside=count,actual_distinct_cap_masks=len(caps),maximum_eligible_mass_exactly='1',receivers=answer)


def geometric_oracle(geom,oracle):
    data=json.loads(geom.read_text()); saved=json.loads(oracle.read_text()); entries=data['caps']+data['witnesses']
    assert all(e['boundary']=='closed' for e in entries) # this production packet uses closed critical witnesses
    bounds=[tuple(map(F,e['bounds'])) for e in entries]; forbidden={tuple(map(F,p)) for p in data['metadata'].get('excluded_points',[])}
    def slack(q):
        cap=sum((F(e['weight']) for e,b in zip(entries,bounds) if e['kind']!='witness' and inside(q,b)),F(0))
        demand=sum((F(e['weight']) for e,b in zip(entries,bounds) if e['kind']=='witness' and inside(q,b)),F(0))
        return cap-demand
    axes=[]
    for k in (0,1):
        ends=sorted({b[j] for b in bounds for j in (2*k,2*k+1)} | {p[k] for p in forbidden})
        reps=[(ends[0]-1,1)]
        for i,a in enumerate(ends):
            reps.append((a,0))
            if i+1<len(ends):reps.append(((a+ends[i+1])/2,1))
        reps.append((ends[-1]+1,1)); axes.append(reps)
    minimum=F(0); negative=0; minimum_by_dim={}; checked=0
    for (x,dx),(y,dy) in product(*axes):
        q=(x,y)
        if q in forbidden:continue
        checked+=1; ss=slack(q)
        if ss<0: negative+=1
        if ss<minimum:minimum=ss
        dim=dx+dy
        if dim not in minimum_by_dim or ss<F(minimum_by_dim[dim]['slack']): minimum_by_dim[dim]=dict(point=fs(q),slack=str(ss))
    assert str(minimum)==saved['minimum_slack']=='-1/3'
    q=tuple(map(F,saved['most_violated_cell']['representative'])); interior=(F(-3,10),F(47,200))
    assert slack(q)==slack(interior)==F(-1,3)
    assert interior[0]<F(-1,25) # outside the exact closed zero-mass collar
    return dict(all_plane_arrangement_atoms_checked=checked,negative_atoms=negative,minimum_slack=str(minimum),minimum_by_dimension=minimum_by_dim,original_reported_point=fs(q),independent_interior_point=fs(interior),both_direct_slacks='-1/3',interior_outside_zero_mass_collar=True,scope='Independent exact full arrangement including unbounded atoms and all forbidden-source coordinates; atom counts differ from producer restricted-union partition.')


def main():
    ap=argparse.ArgumentParser()
    for k in ('model','results','menu','single','augment_model','augment_result','geometry','oracle','output'):ap.add_argument('--'+k.replace('_','-'),type=Path,required=True)
    a=ap.parse_args(); start=time.perf_counter()
    m=json.loads(a.model.read_text()); S=[tuple(F(v,100) for v in p) for p in m['S_positions_hundredths']]; O=[tuple(F(v,100) for v in p) for p in m['outside_positions_hundredths']]
    assert len(S)==25 and len(O)==144 and sum([F(1,25)]*len(S))==1
    centres=[tuple(F(v,100) for v in p) for p in m['centres_hundredths']]; radius=F(m['radius_hundredths'],100)
    assert radius==F(1,4)
    six=json.loads(a.results.read_text())['results'][0]
    rows,metas,counts=load_model(a.model,S,O,centres,radius,True)
    joint=verify_lp(rows,six,len(O)); assert joint['eta']=='-1/100'
    joint['actual_support']=actual_support(S,O,list(map(F,six['primal']['weights'])),centres,radius,False)
    # Check three actual saved single-receiver primals, then D4 transforms
    # giving all six receivers positive feasibility without another solve.
    menu=json.loads(a.menu.read_text())['results']; middle=json.loads(a.single.read_text())['results'][0]
    # Independently reconstruct the saved single-midpoint model and verify its
    # rational dual as well. Only the already source-verified cap row order is
    # reused to interpret the producer's sparse dual row indices.
    single_rows=[row for row,meta in zip(rows,metas) if meta[0]=='cap']
    sx=tuple(F(v,100) for v in middle['centres_hundredths'][0])
    for R in sorted({dist(q,sx) for q in O}):
        hit=[i for i,q in enumerate(O) if dist(q,sx)<=R]
        single_rows.append((hit,[1]*len(hit),4*R*R))
    hit=[i for i,q in enumerate(O) if dist(q,sx)<=radius]
    single_rows.append((hit+[len(O)],[-1]*len(hit)+[1],sum(dist(p,sx)<=radius for p in S)*F(1,25)-4*radius**2))
    single_dual=verify_lp(single_rows,middle,len(O)); assert single_dual['eta']=='7/100'
    bases={tuple(r['centres_hundredths'][0]):r for r in menu+[middle] if len(r['centres_hundredths'])==1 and r['radius_hundredths']==25}
    ids={q:i for i,q in enumerate(O)}; single_feasible=[]
    choices=[((-10,24),'id'),((-10,48),'id'),((-10,0),'flip_y'),((-10,0),'rotate'),((-10,48),'rotate'),((-10,24),'anti')]
    for target,(base,kind) in zip(centres,choices):
        def transform(p):
            x,y=p
            return {'id':(x,y),'flip_y':(x,F(24,25)-y),'rotate':(y,F(24,25)-x),'anti':(F(24,25)-y,F(24,25)-x)}[kind]
        assert transform(tuple(F(v,100) for v in base))==target
        assert set(map(transform,S))==set(S) and set(map(transform,O))==set(O)
        z=[F(0)]*len(O)
        for q,w in zip(O,bases[base]['primal']['weights']):z[ids[transform(q)]]=F(w)
        support=actual_support(S,O,z,[target],radius,True)
        assert support['receivers'][0]['closed_witness_margin']==bases[base]['primal']['eta']
        single_feasible.append(dict(transformed_from=list(base),D4=kind,**support))
    aug=json.loads(a.augment_result.read_text()); AO=[tuple(map(F,p)) for p in aug['primal']['positions']]
    assert AO[:len(O)]==O and AO[len(O):]==[(F(-3,10),F(47,200))]
    ar,am,ac=load_model(a.augment_model,S,AO,centres,radius,False)
    augmented=verify_lp(ar,aug,len(AO)); assert augmented['eta']=='101/10000'
    augmented['actual_support']=actual_support(S,AO,list(map(F,aug['primal']['weights'])),centres,radius,True)
    oracle=geometric_oracle(a.geometry,a.oracle)
    result=dict(status='PASS_INDEPENDENT_R71_EXACT_PRIMAL_DUAL_SOURCE_AND_ORACLE_AUDIT',source_sha256=sha(Path(__file__)),python=sys.version,input_sha256={k:sha(getattr(a,k)) for k in ('model','results','menu','single','augment_model','augment_result','geometry','oracle')},six_model=counts,six_certificate=joint,single_midpoint_exact_certificate=single_dual,six_individual_feasible=single_feasible,one_column_model=ac,one_column_certificate=augmented,continuous_dual=oracle,elapsed_seconds=time.perf_counter()-start,scope='No LP solve or production imports. Independently reconstruct all coordinate rectangle minimum caps, residual ties and row signs; exact primal/dual optimum checks for original six, saved single midpoint and one-column models; six single primals via exact D4 of three saved results; direct positive support and larger open witnesses; independent all-plane arrangement. No area or general C2 bound.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],six_eta=joint['eta'],augmented_eta=augmented['eta'],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
