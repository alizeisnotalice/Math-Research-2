#!/usr/bin/env python3
"""R71 bounded shared-outside LP; exact integer geometry, floating proposals.

Positions and receiver coordinates use hundredths. All rhs use denominator 10000.
Closed captures suffice for the residual non-strict <= condition. A strict full
witness is certified only by an exactly positive rational margin.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, sys, time
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
from scipy.optimize import linprog
from scipy.io import savemat

D=10000
S=[(24*i,24*j) for i in range(5) for j in range(5)]
O=[(24*i,24*j) for i in range(-4,9) for j in range(-4,9) if not (0<=i<=4 and 0<=j<=4)]

@lru_cache(maxsize=1)
def base():
    rows=[]; rhs=[]; labels=[]
    for il in range(-4,9):
      for ih in range(il,9):
       for jl in range(-4,9):
        for jh in range(jl,9):
         hit=[k for k,(x,y) in enumerate(O) if 24*il<=x<=24*ih and 24*jl<=y<=24*jh]
         if not hit: continue
         ns=sum(24*il<=x<=24*ih and 24*jl<=y<=24*jh for x,y in S)
         rows.append(hit);rhs.append(max(D,576*max(ih-il,jh-jl)**2)-400*ns)
         labels.append(["cap",il,ih,jl,jh,ns])
    return rows,rhs,labels

def in_internal_halo(p):
    ds=sorted(max(abs(x-p[0]),abs(y-p[1])) for x,y in S)
    return any(400*(i+1)>4*d*d for i,d in enumerate(ds))

def make_model(centres, radius=49):
    rr,bb,ll=base();rows=list(rr);rhs=list(bb);labels=list(ll);vals=[[1]*len(r) for r in rows]
    for p in centres:
      dist=[max(abs(x-p[0]),abs(y-p[1])) for x,y in O]
      for d in sorted(set(dist)):
        rows.append([k for k,a in enumerate(dist) if a<=d]); vals.append([1]*len(rows[-1]))
        rhs.append(4*d*d);labels.append(["residual",*p,d])
      hit=[k for k,d in enumerate(dist) if d<=radius]
      ns=sum(max(abs(x-p[0]),abs(y-p[1]))<=radius for x,y in S)
      rows.append(hit+[len(O)]);vals.append([-1]*len(hit)+[1])
      rhs.append(400*ns-4*radius*radius);labels.append(["witness",*p,radius,ns])
    ri=[];ci=[];vv=[]
    for i,(rr,vs) in enumerate(zip(rows,vals)):
      ri.extend([i]*len(rr));ci.extend(rr);vv.extend(vs)
    A=sparse.csr_matrix((vv,(ri,ci)),shape=(len(rows),len(O)+1),dtype=float)
    return A,np.array(rhs,dtype=float)/D,rows,vals,rhs,labels

def exact_primal(v,rows,vals,rhs):
    q=[F(float(a)).limit_denominator(10**8) for a in v[:-1]]
    scale=F(1)
    for rr,vs,b in zip(rows,vals,rhs):
      if len(O) in rr: continue
      z=sum(q[j]*a for j,a in zip(rr,vs))
      if z>F(b,D): scale=min(scale,F(b,D)/z)
    q=[a*scale for a in q]
    margins=[]; max_violation=F(0)
    for rr,vs,b in zip(rows,vals,rhs):
      z=sum(q[j]*a for j,a in zip(rr,vs) if j<len(O))
      if len(O) in rr:margins.append(F(b,D)-z)
      else:max_violation=max(max_violation,z-F(b,D))
    eta=min(margins)
    return {"weights":[str(a) for a in q],"eta":str(eta),"eta_float":float(eta),"scale":str(scale),"max_exact_nonwitness_violation":str(max_violation),"positive_positions":sum(a>0 for a in q)}

def exact_dual(res,A,rhs,labels):
    # linprog solves min(-eta): its inequality marginals are nonpositive.
    y=[F(max(0.,float(-a))).limit_denominator(10**7) for a in res.ineqlin.marginals]
    witness=[i for i,l in enumerate(labels) if l[0]=="witness"]
    sw=sum(y[i] for i in witness)
    if sw<=0:return {"status":"no_witness_weight"}
    y=[a/sw for a in y]
    col=[F(0) for _ in O]
    coo=A.tocoo()
    for i,j,a in zip(coo.row,coo.col,coo.data):
      if j<len(O):col[j]+=int(a)*y[i]
    # Every outside singleton is one of the cap rows. Add exact cover deficits.
    added=[]
    for j,c in enumerate(col):
      if c>=0:continue
      x,ypos=O[j]; lab=["cap",x//24,x//24,ypos//24,ypos//24,0]
      i=labels.index(lab); y[i]+=-c;added.append([j,str(-c)])
    bound=sum(a*F(b,D) for a,b in zip(y,rhs))
    return {"status":"exact_feasible_dual", "upper_eta":str(bound),"upper_eta_float":float(bound),"nonzero_rows":[[i,str(a)] for i,a in enumerate(y) if a],"deficit_repairs":added,"sum_witness":"1"}

def solve(centres,radius=49,export=None):
    A,b,rows,vals,rhs,labels=make_model(centres,radius)
    c=np.zeros(len(O)+1);c[-1]=-1
    t=time.monotonic();res=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None)]*len(O)+[(None,None)],method="highs",options={"dual_feasibility_tolerance":1e-9,"primal_feasibility_tolerance":1e-9})
    out={"centres_hundredths":centres,"radius_hundredths":radius,"internal_halo":[in_internal_halo(p) for p in centres],"rows":A.shape[0],"columns":A.shape[1],"nnz":A.nnz,"status":res.status,"message":res.message,"seconds":time.monotonic()-t}
    if res.success:
      out.update(eta_float=float(res.x[-1]),max_float_violation=float(np.max(A@res.x-b)),primal=exact_primal(res.x,rows,vals,rhs),dual=exact_dual(res,A,rhs,labels))
    if export:
      savemat(str(export)+".mat",{"A":A,"b":b,"f":c[:,None],"lb":np.r_[np.zeros(len(O)),-np.inf][:,None],"ub":np.full((len(O)+1,1),np.inf),"outside_positions_hundredths":np.array(O),"S_positions_hundredths":np.array(S),"centres_hundredths":np.array(centres),"radius_hundredths":radius,"scipy_x":res.x[:,None],"scipy_objective":res.fun},do_compression=True)
      Path(str(export)+"-exact.json").write_text(json.dumps({"rhs_denominator":D,"rhs_numerators":rhs,"rows":rows,"values":vals,"labels":labels,"outside_positions_hundredths":O,"S_positions_hundredths":S,"centres_hundredths":centres,"radius_hundredths":radius},indent=2)+"\n")
      d={"schema":"r071-continuous-dual-v1","witnesses":[],"caps":[],"metadata":{"S":[{"position":[str(F(x,100)),str(F(y,100))],"mass":"1/25"} for x,y in S],"source":"joint_outside_lp.py", "finite_grid_only":True,"excluded_points":[[str(F(x,100)),str(F(y,100))] for x,y in S]}}
      for i,a in out["dual"]["nonzero_rows"]:
        l=labels[i]; kind=l[0]
        if kind=="cap":
          _,il,ih,jl,jh,_=l;bounds=[F(24*il,100),F(24*ih,100),F(24*jl,100),F(24*jh,100)]
        else:
          _,x,y,r,*_=l;bounds=[F(x-r,100),F(x+r,100),F(y-r,100),F(y+r,100)]
        entry={"id":str(i),"bounds":[str(z) for z in bounds],"weight":a,"boundary":"closed","rhs":str(F(rhs[i],D)),"kind":kind}
        d["witnesses" if kind=="witness" else "caps"].append(entry)
      Path(str(export)+"-dual-geometry.json").write_text(json.dumps(d,indent=2)+"\n")
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True);ap.add_argument("--export");ap.add_argument("--centres");ap.add_argument("--radius",type=int,default=49);ap.add_argument("--menu",action="store_true");args=ap.parse_args()
    t=time.monotonic()
    if args.centres: requests=[json.loads(args.centres)]
    else:requests=[[[x,y]] for x,y in [(-40,-40),(-30,-30),(-20,-20),(-10,-10),(-40,48),(-30,48),(-20,48),(-10,48),(-40,0),(-30,0),(-20,0),(-10,0),(-40,24),(-20,24)]]
    if args.menu:
      requests=[[[x,y]] for x,y in [(-10,0),(-20,0),(-10,24),(-20,24),(-10,-10),(-20,-20),(-30,0),(-40,0)]]
      results=[solve(ps,radius=r) for ps in requests for r in [15,25,37,49]]
    else:results=[solve(ps,radius=args.radius,export=args.export if len(requests)==1 else None) for ps in requests]
    out={"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform(),"pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"rhs_denominator":D,"S_mass":"1","S_side":"24/25","outside_count":len(O),"seed":None,"scope":"fixed S25, fixed outside grid; fixed witness radius 49/100; all cap rectangles and all residual distance events","results":results,"seconds":time.monotonic()-t}
    Path(args.out).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"output":args.out,"pid":os.getpid(),"seconds":out["seconds"],"models":[{"p":r["centres_hundredths"],"r":r["radius_hundredths"],"eta":r.get("eta_float"),"exact_eta":r.get("primal",{}).get("eta_float"),"exact_dual":r.get("dual",{}).get("upper_eta_float"),"internal_halo":r["internal_halo"]} for r in results]},indent=2))

if __name__=="__main__":main()
