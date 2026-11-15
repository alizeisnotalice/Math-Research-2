#!/usr/bin/env python3
"""Independent R52 saved-bbox mass binding and endpoint-jump union replay.
No candidate enumeration, production derivative code or UnionTree is imported.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import gzip
import json
from math import isqrt
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode=True

def norm(value):
    if isinstance(value,Q):return str(value)
    if isinstance(value,dict):return {str(k):norm(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [norm(v) for v in value]
    return value

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def signature(value):return hashlib.sha256(json.dumps(norm(value),sort_keys=True,separators=(',',':')).encode()).hexdigest()

@lru_cache(None)
def sqrtbox(q,bits):
    scale=1<<bits;k=isqrt((q.numerator<<(2*bits))//q.denominator)
    lo=Q(k,scale);return lo,lo if lo*lo==q else Q(k+1,scale)

def canonical(row):
    a,v,s,u,b=row;a,v,u,b=map(Q,(a,v,u,b));s=int(s)
    if not s:return a,v,0,Q(0),Q(0)
    n,d=isqrt(u.numerator),isqrt(u.denominator)
    if not b and n*n==u.numerator and d*d==u.denominator:return a+s*Q(n,d),v,0,Q(0),Q(0)
    return a,v,s,u,b

def at0(e,bits):
    a,v,s,u,b=e;l,r=sqrtbox(u,bits)
    return (a+s*l,a+s*r) if s>=0 else (a+s*r,a+s*l)

def union_ends(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('source-slab','joint-slab','center-data','mass-velocity','certificate','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();begin=time.process_time()
    def read(path):
        raw=path.read_bytes();return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)
    old,joint,center,mv,cert=[read(x) for x in (a.source_slab,a.joint_slab,a.center_data,a.mass_velocity,a.certificate)]
    assert sha(a.joint_slab)==cert['slab_file_sha256']
    assert cert['source_candidate_code_sha256']==joint['code_sha256']
    assert signature(center)==joint['center_data_canonical_sha256']==old['input_canonical_sha256']
    assert signature(list(map(Q,mv)))==joint['mass_velocity_sha256']
    classes=[(Q(x),int(v)) for x,v in old['x_classes']];ys=list(map(Q,old['y_classes']))
    nx,ny=len(classes),len(ys);ix={v:i+1 for i,v in enumerate(classes)};iy={v:i+1 for i,v in enumerate(ys)}
    pts=[list(map(Q,row)) for row in center['pts']];vel=list(map(Q,mv));z0=Q(joint['z0']);Px,Py=Q(center['Px']),Q(center['Py'])
    P=[[Q(0)]*(ny+1) for _ in range(nx+1)];D=[[Q(0)]*(ny+1) for _ in range(nx+1)]
    sourcechecks=0
    for src in old['source_catalog']:
        label=src['label'];x,y,w=pts[label];i,j=src['period_shift'];v=old['velocities'][label]
        assert (Q(src['x']),Q(src['y']),Q(src['w']),src['v'])==(x+v*z0+i*Px,y+j*Py,w,v)
        s,t=ix[Q(src['x']),v],iy[Q(src['y'])]
        P[s][t]+=w;D[s][t]+=vel[label];sourcechecks+=1
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            for prefix in (P,D):prefix[i][j]+=prefix[i-1][j]+prefix[i][j-1]-prefix[i-1][j-1]
    def extract(prefix,key):
        i,k,j,l=key;return prefix[k+1][l+1]-prefix[i][l+1]-prefix[k+1][j]+prefix[i][j]
    assert {tuple(r['bbox']) for r in old['box_records']}=={tuple(r['bbox']) for r in joint['box_records']}
    lam=Q(joint['lambda']);R=Q(joint['R']);reconstructed=[];nonzero=0;squares=0
    for record,raw in zip(joint['box_records'],joint['boxes']):
        key=record['bbox'];i,k,j,l=key;m,dm=extract(P,key),extract(D,key);q,b=m/(4*lam),dm/(4*lam)
        assert m==Q(record['m0']) and dm==Q(record['mu']) and q==Q(record['q0']) and b==Q(record['q1'])
        assert record['radius_branch']=='root' and q<R*R
        expected=[(classes[k][0],classes[k][1],-1,q,b),(classes[i][0],classes[i][1],1,q,b),
                  (ys[l],0,-1,q,b),(ys[j],0,1,q,b)]
        assert [canonical(e) for e in raw]==[canonical(e) for e in expected]
        reconstructed.append([key,q,b]);nonzero+=bool(b)
        n,d=isqrt(q.numerator),isqrt(q.denominator)
        squares+=bool(b and n*n==q.numerator and d*d==q.denominator)
    bits=cert['bits'];boxes=[tuple(canonical(e) for e in box) for box in joint['boxes']]
    xl,xr,yl,yr=map(Q,joint['domain']);fixed=lambda q:canonical([q,0,0,0,0]);XL,XR,YL,YR=map(fixed,(xl,xr,yl,yr))
    X={XL,XR};Y={YL,YR}
    for l,r,b,t in boxes:X.update((l,r));Y.update((b,t))
    xrows=sorted((at0(e,bits)[0],at0(e,bits)[1],e) for e in X);yrows=sorted((at0(e,bits)[0],at0(e,bits)[1],e) for e in Y)
    assert all(l[1]<r[0] for rows in (xrows,yrows) for l,r in zip(rows,rows[1:]))
    assert signature([e for l,h,e in xrows])==cert['x_endpoint_order_sha256']
    assert signature([e for l,h,e in yrows])==cert['y_endpoint_order_sha256']
    # Independently check final rectangular radius: center gaps dominate full
    # dz drift plus dt mean-value drift, and remain at least half center gap.
    Z,T=Q(cert['dz_radius']),Q(cert['dt_radius']);BT=Q(cert['base_dt_radius']);checked=0
    for rows in (xrows,yrows):
        for left,right in zip(rows,rows[1:]):
            drift=abs(left[2][1]-right[2][1])*Z
            for e in (left[2],right[2]):
                if e[2]:
                    qlo=e[3]-abs(e[4])*BT;assert qlo>0
                    rlo=sqrtbox(qlo,bits)[0];assert rlo>0
                    drift+=abs(e[4])*T/(2*rlo)
            gap=right[0]-left[1];assert 2*drift<=gap;checked+=1
    # Independent telescope: beta=sum v_e*(Y_before(e)-Y_after(e)).
    # Production instead uses adjacent-gap velocity differences times Y_gap.
    rx={e:i for i,(_,_,e) in enumerate(xrows)};ry={e:i for i,(_,_,e) in enumerate(yrows)}
    ev=[[] for _ in xrows];yb=[]
    for i,(l,r,b,t) in enumerate(boxes):
        lo,hi=max(rx[l],rx[XL]),min(rx[r],rx[XR]);bot,top=max(ry[b],ry[YL]),min(ry[t],ry[YR]);assert lo<hi and bot<top
        yb.append((bot,top));ev[lo].append((1,i));ev[hi].append((-1,i))
    coefficients=defaultdict(Q);active=set();jumps=0
    def add_union(weight):
        for lo,hi in union_ends(yb[i] for i in active):
            coefficients[yrows[hi][2]]+=weight;coefficients[yrows[lo][2]]-=weight
    for j,(_,_,e) in enumerate(xrows):
        v=e[1]
        if v:add_union(v);jumps+=1
        for sign,i in ev[j]:
            if sign==1:assert i not in active;active.add(i)
            else:assert i in active;active.remove(i)
        if v:add_union(-v)
    assert not active
    coefficients={e:c for e,c in coefficients.items() if c}
    assert [[e,c] for e,c in sorted(coefficients.items())]==[[canonical(e),Q(c)] for e,c in cert['nonzero_y_endpoint_coefficients']]
    const=Q(0);terms=defaultdict(Q)
    for e,c in coefficients.items():
        const+=c*e[0]
        if e[2]:terms[e[3],e[4]]+=c*e[2]
    terms={k:c for k,c in terms.items() if c}
    assert const==Q(cert['beta_constant'])
    assert [{'q0':u,'q1':b,'coefficient':c} for (u,b),c in sorted(terms.items())]==[{k:Q(v) for k,v in term.items()} for term in cert['beta_root_terms']]
    # Center beta bounds must refine the frozen 80-bit one-axis certificate.
    old_slope_path=a.source_slab.with_name('slope.json');old_slope=json.loads(old_slope_path.read_text())
    bl,bh=map(Q,cert['beta_center_bounds']);assert Q(old_slope['beta_lower'])<=bl<=bh<=Q(old_slope['beta_upper'])
    assert sum(p[2] for p in pts)==Q(cert['mass']) and sum(vel)==Q(cert['mass_t_derivative'])
    # A raw crossing can be entirely masked: A=[-2,-s]x[0,1],
    # B=[s,2]x[0,1], C=[-1,1]x[0,1], |s|<1/2. The raw facing
    # endpoints swap at0 while the union is identically [-2,2]x[0,1].
    for s in map(Q,['-1/4','0','1/4']):assert union_ends([(-2,-s),(s,2),(-1,1)])==[(-2,2)]
    result={'status':'passed','scope':'independent source-prefix binding for every supplied active bbox plus direct endpoint-jump union replay; complete rejected-bbox enumeration remains separate premise',
      'source_identity_checks':sourcechecks,'all_active_bbox_mass_and_derivative_bindings':len(reconstructed),
      'nonzero_bbox_mass_derivatives':nonzero,'moving_rational_square_roots':squares,'bbox_q0_q1_sha256':signature(reconstructed),
      'final_radius_consecutive_gap_checks':checked,'moving_endpoint_jump_checks':jumps,
      'root_sum_matches_gap_producer_exactly':True,'center_beta_refines_frozen_one_axis_certificate':True,
      'masked_raw_wall_example':'A=[-2,-s]x[0,1],B=[s,2]x[0,1],C=[-1,1]x[0,1]; for |s|<1/2 union=[-2,2]x[0,1], although raw facing endpoints -s and s swap',
      'source_sha256':{k:sha(v) for k,v in vars(a).items() if k!='output'},'old_slope_sha256':sha(old_slope_path),
      'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-begin,
      'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,*sys.argv]}
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(norm(result),indent=2)+'\n');print(json.dumps(norm(result),indent=2))
if __name__=='__main__':main()
