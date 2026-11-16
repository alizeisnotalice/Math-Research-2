#!/usr/bin/env python3
"""One exact-column R71 support augmentation, rebuilding all finite caps.

All positions use units 1/200, rhs denominator 40000. Optional oracle-cell
representatives are read as exact fractions, allowing further bounded steps.
"""
import argparse,hashlib,json,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
from scipy.optimize import linprog
from scipy.io import savemat

def run(extra,stem):
    t=time.monotonic();den=F(40000)
    S=[(F(48*i),F(48*j)) for i in range(5) for j in range(5)]
    O=[(F(48*i),F(48*j)) for i in range(-4,9) for j in range(-4,9) if not(0<=i<=4 and 0<=j<=4)]
    for p in extra:
      q=tuple(F(a)*200 for a in p)
      if q not in S+O:O.append(q)
    n=len(O);xs=sorted(set(x for x,y in S+O));ys=sorted(set(y for x,y in S+O))
    caps={};rectangle_count=0
    for ai,xl in enumerate(xs):
     for xh in xs[ai:]:
      for bi,yl in enumerate(ys):
       for yh in ys[bi:]:
        hit=tuple(k for k,(x,y) in enumerate(O) if xl<=x<=xh and yl<=y<=yh)
        if not hit:continue
        ns=sum(xl<=x<=xh and yl<=y<=yh for x,y in S)
        b=max(den,max(xh-xl,yh-yl)**2)-1600*ns
        rectangle_count+=1
        if hit not in caps or b<caps[hit][0]:caps[hit]=(b,(xl,xh,yl,yh))
    rows=[];vals=[];rhs=[];labels=[]
    for hit,(b,box) in caps.items():rows.append(list(hit));vals.append([1]*len(hit));rhs.append(b);labels.append(["cap",*box])
    # Four genuine side-one zero-cap squares, each containing S. They are
    # redundant on the finite support, but globally valid dual geometry.
    for xl,xh in [(F(-8),F(192)),(F(0),F(200))]:
     for yl,yh in [(F(-8),F(192)),(F(0),F(200))]:
      hit=[k for k,(x,y) in enumerate(O) if xl<=x<=xh and yl<=y<=yh]
      rows.append(hit);vals.append([1]*len(hit));rhs.append(F(0));labels.append(["collar",xl,xh,yl,yh])
    centres=[(-20,48),(-20,96),(-20,192),(0,212),(96,212),(144,212)]
    for p in centres:
      ds=[max(abs(x-p[0]),abs(y-p[1])) for x,y in O]
      for d in sorted(set(ds)):
        hit=[k for k,a in enumerate(ds) if a<=d]
        rows.append(hit);vals.append([1]*len(hit));rhs.append(4*d*d);labels.append(["residual",F(p[0]-d),F(p[0]+d),F(p[1]-d),F(p[1]+d)])
      r=50;hit=[k for k,a in enumerate(ds) if a<=r];ns=sum(max(abs(x-p[0]),abs(y-p[1]))<=r for x,y in S)
      rows.append(hit+[n]);vals.append([-1]*len(hit)+[1]);rhs.append(F(1600*ns-4*r*r));labels.append(["witness",F(p[0]-r),F(p[0]+r),F(p[1]-r),F(p[1]+r)])
    ri=[];ci=[];vv=[]
    for i,(rr,vs) in enumerate(zip(rows,vals)):ri.extend([i]*len(rr));ci.extend(rr);vv.extend(vs)
    A=sparse.csr_matrix((vv,(ri,ci)),shape=(len(rows),n+1),dtype=float);b=np.array([float(a/den) for a in rhs]);c=np.zeros(n+1);c[-1]=-1
    res=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None)]*n+[(None,None)],method="highs",options={"dual_feasibility_tolerance":1e-9,"primal_feasibility_tolerance":1e-9})
    assert res.success,res.message
    q=[F(float(a)).limit_denominator(10**9) for a in res.x[:-1]];scale=F(1)
    for rr,vs,bb in zip(rows,vals,rhs):
      if n in rr:continue
      z=sum(q[j]*a for j,a in zip(rr,vs))
      if z>bb/den:scale=min(scale,(bb/den)/z)
    q=[a*scale for a in q];eta=min(bb/den+sum(q[j] for j in rr if j<n) for rr,bb in zip(rows,rhs) if n in rr)
    maxvio=max([F(0)]+[sum(q[j]*a for j,a in zip(rr,vs))-bb/den for rr,vs,bb in zip(rows,vals,rhs) if n not in rr])
    assert maxvio==0
    y=[F(max(0.,float(-a))).limit_denominator(10**9) for a in res.ineqlin.marginals];sw=sum(y[i] for i,l in enumerate(labels) if l[0]=="witness");assert sw>0;y=[a/sw for a in y]
    col=[F(0)]*n
    for i,a in enumerate(y):
      if a:
       for j,v in zip(rows[i],vals[i]):
        if j<n:col[j]+=v*a
    repairs=[]
    for j,a in enumerate(col):
      if a<0:
        i=next(i for i,rr in enumerate(rows) if rr==[j] and vals[i]==[1]);y[i]-=a;repairs.append([j,str(-a)])
    upper=sum(a*bb/den for a,bb in zip(y,rhs))
    geom={"schema":"r071-continuous-dual-v1","witnesses":[],"caps":[],"metadata":{"S":[{"position":[str(x/200),str(y/200)],"mass":"1/25"} for x,y in S],"extra_positions":extra,"finite_grid_only":True}}
    for i,a in enumerate(y):
      if not a:continue
      l=labels[i];entry={"id":str(i),"bounds":[str(v/200) for v in l[1:]],"boundary":"closed","kind":"cap" if l[0]=="collar" else l[0],"weight":str(a),"rhs":str(rhs[i]/den)}
      geom["witnesses" if l[0]=="witness" else "caps"].append(entry)
    out={"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform(),"pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"outside_count":n,"extra_positions":extra,"centres":[[str(F(x,200)),str(F(y,200))] for x,y in centres],"witness_radius":"1/4","full_rectangle_enumeration":rectangle_count,"deduplicated_cap_count":len(caps),"collar_rows":4,"rows":len(rows),"nnz":A.nnz,"eta_float":res.x[-1],"primal":{"eta":str(eta),"eta_float":float(eta),"scale":str(scale),"weights":[str(a) for a in q],"positions":[[str(x/200),str(y/200)] for x,y in O],"max_exact_nonwitness_violation":str(maxvio)},"dual":{"upper_eta":str(upper),"upper_eta_float":float(upper),"repairs":repairs,"nonzero_rows":[[i,str(a)] for i,a in enumerate(y) if a]},"seconds":time.monotonic()-t,"evidence_scope":"finite fixed support LP only; all bbox caps and complete residual radius events rebuilt; exact rational primal/dual certificates"}
    Path(stem+"-result.json").write_text(json.dumps(out,indent=2)+"\n");Path(stem+"-dual-geometry.json").write_text(json.dumps(geom,indent=2)+"\n")
    Path(stem+"-model.json").write_text(json.dumps({"rhs":[str(a/den) for a in rhs],"rows":rows,"values":vals,"labels":[[l[0]]+[str(a/200) for a in l[1:]] for l in labels]},separators=(",",":"))+"\n")
    savemat(stem+".mat",{"A":A,"b":b[:,None],"f":c[:,None],"lb":np.r_[np.zeros(n),-np.inf][:,None],"ub":np.full((n+1,1),np.inf),"scipy_x":res.x[:,None]},do_compression=True)
    print(json.dumps({k:out[k] for k in ["pid","outside_count","extra_positions","rows","eta_float","seconds"]}|{"exact_eta":str(eta),"exact_upper":str(upper)},indent=2))

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--extra",default='[["-3/10","47/200"]]');ap.add_argument("--stem",required=True);a=ap.parse_args();run(json.loads(a.extra),a.stem)
