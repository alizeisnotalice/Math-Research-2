#!/usr/bin/env python3
"""Independent R54 left transform/bootstrap/full-bbox single-root certificate."""
import argparse, hashlib, json, math, os, sys, time
from fractions import Fraction as F
from pathlib import Path

def digest(p):
    b=p.read_bytes(); return dict(path=str(p),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b))
def ceilq(x): return -(-x.numerator//x.denominator)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--source',type=Path,required=True); ap.add_argument('--transformed',type=Path,required=True); ap.add_argument('--axis',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    assert os.path.ismount('/Volumes/data'); start=time.perf_counter()
    src=json.loads(a.source.read_text()); tr=json.loads(a.transformed.read_text()); ax=json.loads(a.axis.read_text())
    pts=[tuple(map(F,p)) for p in src['pts']]; tp=[tuple(map(F,p)) for p in tr['pts']]
    px,py,lam=(F(src[k]) for k in ('Px','Py','lam')); assert lam==1
    assert all(src[k]==tr[k] for k in ('Px','Py','lam'))
    X=max(x for x,y,m in pts); xmin=min(x for x,y,m in pts)
    assert xmin==0 and X<px and len(pts)==len(tp)==204
    assert tp==[(X-x,y,m) for x,y,m in pts]
    assert ax['input_sha256']==digest(a.transformed)['sha256']
    mass=sum(m for x,y,m in pts); rho=mass/(px*py)
    assert rho*(1+px/20)*(1+py/20)<1
    assert rho*(1+px/15)*(1+py/30)<1
    long=rho*F(3,2)*(1+py/(6*px)); assert long<1 and 6*px>30
    R=F(10); bs=[]
    for step in ax['bootstrap']:
        assert F(step['R'])==R
        nc,nr=ceilq(R/px),ceilq(2*R/py)+1
        row=[(i,k,x+i*px,y,m) for i in range(1-nc,1) for k,(x,y,m) in enumerate(tp) if x+i*px>X-R]
        B=nr*sum(m for i,k,x,y,m in row)
        assert [[i,k] for i,k,x,y,m in row]==step['row_labels']
        assert nc==step['columns'] and nr==step['rows'] and B==F(step['mass_upper'])
        nextR=F(step['next_strict_radius']); assert 4*nextR**2>B
        bs.append(dict(R=str(R),columns=nc,rows=nr,mass_upper=str(B),next_R=str(nextR),positive_square_slack=str(4*nextR**2-B)))
        if nextR>=R: break
        R=nextR
    assert R==F(ax['final_owner_radius']) and nc==1 and nr==3
    atoms={}
    for i,k,x,y,m in row:
        for j in range(nr): atoms[x,y+j*py]=atoms.get((x,y+j*py),F(0))+m
    xs=sorted({x for x,y in atoms}); ys=sorted({y for x,y in atoms}); nx,ny=len(xs),len(ys)
    W=math.lcm(*(m.denominator for m in atoms.values())); xi={x:i+1 for i,x in enumerate(xs)}; yi={y:j+1 for j,y in enumerate(ys)}
    pre=[[0]*(ny+1) for _ in range(nx+1)]
    for (x,y),m in atoms.items():pre[xi[x]][yi[y]]+=int(m*W)
    for i in range(1,nx+1):
        z=0
        for j in range(1,ny+1): z+=pre[i][j]; pre[i][j]=pre[i-1][j]+z
    m0=max(m for (x,y),m in atoms.items() if x==X); q=m0/4
    counts=dict(all_boxes=0,zero=0,empty=0,active=0,p_le_q=0,L_le_zero=0,squared=0,strict=0,equality=0); ties=[]; min_slack=None
    for i,u in enumerate(xs):
      for k,v in enumerate(xs[i:],i):
       for j,c in enumerate(ys):
        for l,d in enumerate(ys[j:],j):
         counts['all_boxes']+=1; M=pre[k+1][l+1]-pre[i][l+1]-pre[k+1][j]+pre[i][j]
         if M==0: counts['zero']+=1; continue
         if max((v-u)**2,(d-c)**2)>=F(M,W): counts['empty']+=1; continue
         counts['active']+=1; p=F(M,4*W); g=X-u; L=p-q-g*g
         if p<=q: counts['p_le_q']+=1; equal=p==q and g==0
         elif L<=0: counts['L_le_zero']+=1; equal=False
         else:
          counts['squared']+=1; margin=4*g*g*q-L*L; assert margin>=0; equal=margin==0
         slack=q+g*g-p
         if slack>0:min_slack=slack if min_slack is None else min(min_slack,slack)
         if equal:counts['equality']+=1;ties.append(dict(bbox=list(map(str,(u,v,c,d))),mass=str(F(M,W))))
         else:counts['strict']+=1
    assert counts['all_boxes']==nx*(nx+1)//2*ny*(ny+1)//2
    assert ties and all(t['bbox'][0]==t['bbox'][1]==str(X) and F(t['mass'])==m0 for t in ties)
    lo,hi=F(ax['c_edge_lower']),F(ax['c_edge_upper']); assert 0<lo-X<hi-X and (lo-X)**2<=q<(hi-X)**2
    assert math.isqrt(q.numerator)**2!=q.numerator or math.isqrt(q.denominator)**2!=q.denominator
    out=dict(status='PASS',python=sys.version,arithmetic='Fraction/full coordinate bbox comparison; no numerical tolerance or random seed',inputs=[digest(a.source),digest(a.transformed),digest(a.axis)],script=digest(Path(__file__)),point_count=len(tp),X=str(X),X_less_Px_margin=str(px-X),axis_catalog=dict(atoms=len(atoms),nx=nx,ny=ny),bootstrap=bs,counts=counts,equality_witnesses=ties,singleton_mass=str(m0),q=str(q),transformed_true_right_axis='X+sqrt(q)',original_true_left_axis='-sqrt(q)',axis_irrational=True,operations_interval_contains_exact_root=True,positive_margins=dict(six_Px_minus_30=str(6*px-30),seven_Px_minus_30=str(7*px-30),long_tail_one_minus_bound=str(1-long)),long_tail_bound=str(long),elapsed_seconds=time.perf_counter()-start,scope='New R54 transformed input; exact left-axis single-root proof and finite-N>=7 locality margins; no new seam area execution')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({k:out[k] for k in ('status','axis_catalog','counts','singleton_mass','q','elapsed_seconds')},indent=2))
if __name__=='__main__':main()
