"""Floating local optimization of one placement of nine tested sets."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.spatial import ConvexHull


def unit(t):
    return np.column_stack((np.cos(t), np.sin(t)))


def shape_boundaries(arc_count, circle_count, smooth_count):
    out = [unit(np.arange(circle_count)*2*np.pi/circle_count)*.5]
    for n in (3,5,7,9):
        radius=1/(2*np.cos(np.pi/(2*n)))
        arcs=[]
        for k in range(n):
            phi=2*np.pi*k/n
            angles=np.linspace(phi+np.pi-np.pi/(2*n),phi+np.pi+np.pi/(2*n),arc_count+1)
            arcs.append(radius*unit(np.asarray([phi]))+unit(angles))
        out.append(np.vstack(arcs))
    for n in (11,19):
        out.append(unit(np.arange(n)*2*np.pi/n)/(2*np.cos(np.pi/(2*n))))
    theta=np.arange(smooth_count)*2*np.pi/smooth_count
    for coefficients in [{3:1/80,5:-19/1200},{5:1/200,15:3/1792}]:
        h=.5+sum(a*np.cos(k*theta) for k,a in coefficients.items())
        hp=-sum(a*k*np.sin(k*theta) for k,a in coefficients.items())
        u=unit(theta)
        out.append(h[:,None]*u+hp[:,None]*np.column_stack((-u[:,1],u[:,0])))
    return out


class SampleHull:
    def __init__(self,arc,circle,smooth):
        self.boundaries=shape_boundaries(arc,circle,smooth)
        self.owners=np.concatenate([np.full(len(q),i) for i,q in enumerate(self.boundaries)])
        self.calls=0

    def points(self,x):
        points=self.boundaries[:1] + [self.boundaries[1]+x[:2]]
        d=2
        for base in self.boundaries[2:]:
            r=x[d]
            co,si=np.cos(r),np.sin(r)
            rot=base@np.array([[co,si],[-si,co]])
            points.append(rot+x[d+1:d+3])
            d+=3
        return np.vstack(points)

    def area_gradient(self,x):
        self.calls+=1
        p=self.points(x)
        hull=ConvexHull(p)
        verts=hull.vertices
        q=p[verts]
        owners=self.owners[verts]
        dif=np.roll(q,-1,axis=0)-np.roll(q,1,axis=0)
        g=np.column_stack((dif[:,1],-dif[:,0]))*.5
        y=np.bincount(owners,weights=g[:,1],minlength=9)
        xx=np.bincount(owners,weights=g[:,0],minlength=9)
        angles=np.zeros(9)
        for j in range(2,9):
            d=2+3*(j-2)
            mask=owners==j
            local=q[mask]-x[d+1:d+3]
            angles[j]=np.sum(g[mask,1]*local[:,0]-g[mask,0]*local[:,1])
        return hull.volume,np.r_[xx[1],y[1],np.stack((angles[2:],xx[2:],y[2:]),axis=1).ravel()]

    def area(self,x):
        self.calls+=1
        return ConvexHull(self.points(x)).volume


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    raw=args.source.read_bytes()
    src=json.loads(raw)
    old=np.array(list(map(float,src["parameter_decimals_interpreted_exactly"])))
    coarse=SampleHull(20,192,512)
    fine=SampleHull(96,768,2048)
    hi=SampleHull(384,4096,8192)
    rng=np.random.default_rng(20261118)
    candidates=[]
    period=2*np.pi/5

    def record(phase,x):
        params=[format(float(v),".17g") for v in x]
        val=hi.area(x)
        candidates.append({"phase":phase,"sample_area":val,"params":params})
        print(phase,format(val,".13g"),flush=True)
        args.output.write_text(json.dumps({
            "source_sha256":sha256(raw).hexdigest(),
            "scipy":scipy.__version__,
            "numpy":np.__version__,
            "seed":20261118,
            "sample_counts":{ "coarse":[20,192,512], "fine":[96,768,2048], "evaluation":[384,4096,8192]},
            "calls":[coarse.calls,fine.calls,hi.calls],
            "records":candidates,
            "best_fixed_eight":min((v for v in candidates if v["phase"].startswith("fixed")),key=lambda v:v["sample_area"]),
            "best":min(candidates,key=lambda v:v["sample_area"]),
        },indent=2)+"\n")

    fixed=[]
    for j in range(12):
        z=np.array([j*period/12,0.,0.])
        simplex=np.vstack((z,z+np.diag([period/24,.003,.003])))
        result=minimize(lambda z:coarse.area(np.r_[old,z]),z,method="Nelder-Mead",
                        options={"initial_simplex":simplex,"maxiter":350,"fatol":1e-12,"xatol":1e-8})
        v=np.r_[old,result.x]
        fixed.append((fine.area(v),j,v))
    for _,j,v in sorted(fixed)[:4]:
        opts=minimize(lambda z:fine.area_gradient(np.r_[old,z])[0],v[-3:],
                      method="Nelder-Mead",options={"maxiter":250,"fatol":1e-13,"xatol":1e-8,
                      "initial_simplex":np.vstack((v[-3:],v[-3:]+np.diag([.003,.0002,.0002])))})
        record(f"fixed-{j}",np.r_[old,opts.x])
    scale=np.array([.018,.018]+7*[.1,.018,.018])
    for _,j,seed in sorted(fixed)[:4]:
        for k in range(3):
            x0=seed+(rng.normal(0,0.025,len(seed))*scale if k else 0)
            a=minimize(lambda z:(lambda vp:(vp[0],vp[1]*scale))(coarse.area_gradient(x0+scale*z)),
                       np.zeros(len(x0)),jac=True,method="L-BFGS-B",bounds=[(-1.5,1.5)]*len(x0),
                       options={"maxiter":550,"gtol":1e-10,"ftol":1e-14,"maxls":45})
            x1=x0+scale*a.x
            b=minimize(lambda z:(lambda vp:(vp[0],vp[1]*scale))(fine.area_gradient(x1+scale*z)),
                       np.zeros(len(x1)),jac=True,method="L-BFGS-B",bounds=[(-.3,.3)]*len(x1),
                       options={"maxiter":420,"gtol":1e-11,"ftol":1e-14,"maxls":45})
            x2=x1+scale*b.x
            record(f"joint-{j}-{k}",x2)
    best=min(candidates,key=lambda v:v["sample_area"])
    start=np.asarray(list(map(float,best["params"])))
    for i in range(2):
        r=minimize(lambda z:fine.area(start+scale*z),np.zeros(len(start)),method="Powell",
                   options={"maxfev":9000,"xtol":1e-8,"ftol":1e-12})
        start=start+scale*r.x
        record(f"refine-{i}",start)


if __name__=="__main__":
    main()
