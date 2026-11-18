"""Suggest fixed-label point sets by convex cone separation, then certify with Arb."""
from pathlib import Path
from itertools import combinations
from argparse import ArgumentParser
from hashlib import sha256
import platform
import json


import numpy as np
import scipy
from scipy.optimize import linprog

from certify_four_body_contact_stationary import BOX, geometry, interval, q
from certify_multi_host_points import LIMITS, SEED_DENOM, lower, points, replay, two_hulls

def poses(raw):
    data=json.loads(raw)['zero_gain_five']
    return (
        np.array([float(x.mid()) for x in geometry(*(interval(*row) for row in BOX))[0]]),
        np.array(list(map(float,data['params']))),
    )


pi=np.pi
u=lambda a:np.stack((np.cos(a),np.sin(a)),axis=-1)
m=lambda a:np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
opts={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9}
den=10**12

def body_support(n, phi, theta):
    turn=2*pi/n
    delta=(theta-phi+turn/2)%turn-turn/2
    v=np.abs(delta)<turn/4
    rr=1/(2*np.cos(pi/(2*n)))
    return np.where(v,rr*np.cos(delta),1-rr*np.cos(turn/2-np.abs(delta)))

def host_support(theta, p):
    normal=u(theta)
    h=np.full(len(theta),.5)
    h=np.maximum(h,body_support(3,0,theta)+normal@p[:2])
    for n,j in zip((5,7,9),(2,5,8)):
        if j<len(p):
            h=np.maximum(h,body_support(n,p[j],theta)+normal@p[j+1:j+3])
    return h

class Sampler:
    def __init__(self, host_poses, cells, dirs):
        self.cells=cells
        self.half=pi/cells
        self.angles=(np.arange(cells)+.5)*2*pi/cells
        base=np.arange(dirs)*2*pi/dirs
        self.normals=[]
        self.host_table=[]
        for p in host_poses:
            peaks=[0.512504445326185, 1.577473015380595, 5.763953958980647,
                   .513440833040228,1.5782924138567,2.7783,3.15]
            aa=np.unique(np.r_[base,peaks,np.array(peaks)+pi])
            U=u(aa)
            self.normals.append(aa)
            self.host_table.append((U,host_support(aa,p),np.c_[U,-np.ones(len(U))]))
    def angles_rows(self, points, which=None, fixed_angles=None):
        rows=[]
        for H,(U,host,coef) in enumerate(self.host_table):
            if which is not None and H!=which:continue
            for sign in (1,-1):
                for i,a in enumerate(self.angles if fixed_angles is None else fixed_angles):
                    M=m(a)@np.diag([1,sign])
                    L=U@M
                    mat=L@points.T
                    at=np.argmax(mat,axis=1)
                    maxima=mat[np.arange(len(U)),at]
                    sol=linprog([0,0,1],A_ub=coef,b_ub=host-maxima,bounds=[(None,None)]*3,
                                method='highs',options=opts)
                    assert sol.success
                    w=-sol.ineqlin.marginals
                    jj=np.flatnonzero(w>1e-9)
                    if len(jj)!=3:
                        # Opposite translation facets can become the only optimal pair
                        # for an already fitting point set. A third adjacent facet is
                        # selected from the active face by a tiny objective perturbation.
                        sol=linprog([1.7e-8,2.1e-8,1],A_ub=coef,b_ub=host-maxima,
                                    bounds=[(None,None)]*3,method='highs',options=opts)
                        w=-sol.ineqlin.marginals
                        jj=np.flatnonzero(w>1e-9)
                    if len(jj)!=3:continue
                    W=np.array([np.linalg.det(U[jj[[1,2]]]),np.linalg.det(U[jj[[2,0]]]),np.linalg.det(U[jj[[0,1]]])])
                    if (W<0).all():W*=-1
                    if not (W>1e-10).all():continue
                    W/=W.sum()
                    rows.append((H,sign,i,float(a),jj,at[jj],W,float(sol.fun),sol.x[:2].tolist()))
        return rows

def make_regular(n):
    p=u(np.arange(n)*2*pi/n)/(2*np.cos(pi/(2*n)))
    return p-p[0]

def conic_outer(P,sampler,rows,epsilon=2e-9):
    k=len(P)
    constraints=[]
    rhs=[]
    for H,sign,i,a,jj,at,W,_,t in rows:
        U,h,_=sampler.host_table[H]
        for b in (-1,1):
            M=m(a+b*sampler.half)@np.diag([1,sign])
            v=np.zeros((k,2))
            np.add.at(v,at,W[:,None]*(U[jj]@M))
            constraints.append(np.r_[-v.ravel(),1])
            rhs.append(-W@h[jj])
    pair=list(combinations(range(k),2))
    def cut(a,b,u):
        v=np.zeros(2*k+1)
        v[2*a:2*a+2]=u
        v[2*b:2*b+2]=-u
        constraints.append(v)
        rhs.append(1.0)
    for a,b in pair:
        for d in (0.,pi/2,pi,1.5*pi):
            cut(a,b,u(d))
    for _ in range(70):
        objective=np.r_[np.zeros(2*k),-1]
        lp=linprog(objective,A_ub=np.array(constraints),b_ub=np.array(rhs),bounds=[(0,0)]*2+[(-1,1)]*(2*k-2)+[(None,None)],
                   method='highs',options=opts)
        assert lp.success,lp.message
        p=lp.x[:-1].reshape(-1,2)
        found=0
        maxdist=0
        for a,b in pair:
            d=p[a]-p[b];dist=np.linalg.norm(d);maxdist=max(maxdist,dist)
            if dist>1+epsilon:
                cut(a,b,d/dist)
                found+=1
        if not found: break
    return p / (max(1,maxdist)*(1+1e-11)),float(lp.x[-1]),maxdist,_

def obstruction(rows):
    return {f'{h}/{s}':min((r[7] for r in rows if r[0]==h and r[1]==s),default=-1)
            for h in (0,1) for s in (1,-1)}

def seed(v):
    x,y=v
    s=1 if x>=0 else -1
    return [s,round(SEED_DENOM*s*y/(1+s*x))]


def choose(tab, P, sg, a):
    U,b,F=tab
    M=m(a)@np.diag([1,sg])
    dots=(U@M)@P.T
    h=dots.max(axis=1)
    lp=linprog([0,0,1],A_ub=F,b_ub=b-h,bounds=[(None,None)]*3,method='highs',options=opts)
    assert lp.success
    idx=np.flatnonzero(-lp.ineqlin.marginals>1e-9)
    if len(idx)!=3:
        lp=linprog([1.7e-9,2.1e-9,1],A_ub=F,b_ub=b-h,bounds=[(None,None)]*3,
                   method='highs',options=opts)
        idx=np.flatnonzero(-lp.ineqlin.marginals>1e-9)
        if len(idx)!=3:return [],lp.fun
    return [[*seed(U[j]),int(np.argmax(dots[j]))] for j in idx],lp.fun


def suggest(host, tab, P, exact, s, target):
    final=[]
    stack=[(2,i) for i in range(3,-1,-1)]
    while stack:
        d,j=stack.pop()
        a=(2*j+1)*float(pi)/2**d
        facets,value=choose(tab,P,s,a)
        if facets:
            try:
                lb=lower(host,exact,s,d,j,facets)
            except AssertionError:
                pass
            else:
                if lb>q(target):
                    final.append([d,j,facets])
                    continue
        if d>=19:
            print('failed',float(a),value,flush=True)
            return None
        stack.extend([(d+1,2*j+1),(d+1,2*j)])
    return final

def main():
    ap=ArgumentParser(description="Generate the eleven-point SOCP candidate and strict two-host certificates.")
    ap.add_argument('competing_fits',type=Path)
    ap.add_argument('certificate',type=Path)
    args=ap.parse_args()
    raw=args.competing_fits.read_bytes()
    host_poses=poses(raw)
    candidate=make_regular(11)
    best=np.array(candidate)
    best_margin=-1.
    history=[]
    sparse=Sampler(host_poses,180,220)
    for attempt in range(3):
        active=sparse.angles_rows(candidate)
        stats=obstruction(active)
        margin=min(stats.values())
        if margin>best_margin:
            best_margin=margin
            best=candidate.copy()
        summary={'attempt':attempt,'centre_margins':stats}
        if attempt<2:
            candidate,eta,distance,cuts=conic_outer(candidate,sparse,active)
            summary.update(conic_objective=eta,max_pair_before_rescale=distance,separation_iterations=cuts)
        history.append(summary)
    fine=Sampler(host_poses,720,440)
    active=fine.angles_rows(best)
    P,eta,distance,cuts=conic_outer(best,fine,active)
    p=(P/(1+1e-10)*den).round().astype('int64')
    cert=dict(point_denominator=den,points=p.tolist(),competing_fits_input_sha256=sha256(raw).hexdigest(),
              padding_thresholds={k:str(v) for k,v in LIMITS.items()},proofs={})
    exact,_=points(cert)
    hosts=two_hulls(raw)
    verification=Sampler(host_poses,720,720)
    for j,(host_name,host) in enumerate(hosts.items()):
        cert['proofs'][host_name]={}
        for hand in (1,-1):
            rows=suggest(host,verification.host_table[j],p/den,exact,hand,LIMITS[host_name])
            if rows is None:
                raise ValueError('The numerical optimum did not give a strict witness with these thresholds.')
            cert['proofs'][host_name][str(hand)]=rows
    receipt=replay(raw,cert)
    encoded=(json.dumps(cert,separators=(',',':'))+'\n').encode()
    args.certificate.write_bytes(encoded)
    receipt['search']={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                      'coarse':history,'final_conic_endpoint_objective':eta,
                      'final_max_pair_before_rescale':distance,'separation_iterations':cuts}
    receipt['certificate_sha256']=sha256(encoded).hexdigest()
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
