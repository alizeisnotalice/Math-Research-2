"""Compute forced support offsets from exact width-one, root-isolated curvature profiles."""
from argparse import ArgumentParser
from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,arb_mat,ctx


def rational(s):
    q=Fraction(s)
    return arb(q.numerator)/q.denominator

def step(z,f,t):
    s,c=t.sin_cos()
    y,v=z
    return f+(y-f)*c+v*s, -(y-f)*s+v*c


class Profile:
    def __init__(self,spec):
        self.grid=spec['grid']
        self.bounds=[0]+[k for k,f in spec['pieces']]
        assert self.bounds==sorted(set(self.bounds)) and self.bounds[-1]==self.grid and self.grid%3==0
        self.knots=[k*arb.pi()/self.grid for k in self.bounds]
        self.base=[rational(f) if f is not None else None for k,f in spec['pieces']]
        assert sum(t is None for t in self.base)==3
        def coefficients(x):
            var=iter(x[:3])
            return [next(var) if v is None else v for v in self.base]
        def evaluate(x):
            fs=coefficients(x);states=[(arb(0),x[3])]
            for j,f in enumerate(fs):
                states.append(step(states[-1],f,self.knots[j+1]-self.knots[j]))
            return fs,states
        def constraints(x):
            self.fs,self.states=evaluate(x)
            return [self.at(1,3)[0],self.at(2,3)[0],self.states[-1][0],self.states[-1][1]+x[3]]
        z=constraints([arb(0)]*4)
        e=[constraints([arb(i==j) for i in range(4)])for j in range(4)]
        mat=arb_mat([[e[j][i]-z[i] for j in range(4)]for i in range(4)])
        self.det=mat.det()
        assert abs(self.det)>rational('1e-12')
        result=mat.solve(arb_mat([[-u]for u in z]))
        self.x=[result[i,0] for i in range(4)]
        assert all(v.contains(0) for v in constraints(self.x))
        assert abs(self.x[3])<rational('.2')
        assert all(v>=-arb(1)/2 and v<=arb(1)/2 for v in self.fs)
    def at(self,k,den):
        # Input angle pi*k/den. The antipodal extension is exact.
        k%=2*den
        sgn=1 if k<den else -1
        k%=den
        seg=bisect_right(self.bounds,(k*self.grid)//den)-1
        y,v=step(self.states[seg],self.fs[seg],arb.pi()*k/den-self.knots[seg])
        return sgn*y,sgn*v
    def root_clearance(self):
        m=2000
        dx=arb.pi()/(3*m)
        def F(k,der):
            return sum((-1)**j*self.at(k+m*j,3*m)[int(der)] for j in range(3))
        margin=-(max(F(k,False).upper() for k in range(1,m))+6*dx**2/8)
        ends=min((-F(0,True)-6*dx).lower(),(F(m,True)-6*dx).lower())
        assert margin>rational('1e-7') and ends>rational('1e-5')
        return str(margin),str(ends)
    def support(self,u):
        pi=arb.pi()
        theta=arb.atan2(u[1],u[0])
        if theta<0:theta+=2*pi
        sign=1
        if theta>pi:theta-=pi;sign=-1
        assert 0<theta<pi
        out=None
        for j,f in enumerate(self.fs):
            if theta<self.knots[j] or theta>self.knots[j+1]:continue
            val=arb(1)/2+sign*step(self.states[j],f,theta-self.knots[j])[0]
            out=val if out is None else out.union(val)
        assert out is not None
        return out


def run(comp,source):
    sys.path.insert(0,str(comp.resolve()))
    from certify_ah_tapered_minkowski_closure import geometry
    prepared,_=geometry(comp)
    ctx.prec=192
    allrows={a:(u,h)for a,u,h in prepared['active']}
    cap_rows=['K-C3','deep','I-E3','I-E3','I-E3','I-E3','I-E3','I-E3',None,'K-C3','K-C3','K-C3']
    poses=[(k,s)for k in range(0,360,60)for s in (-1,1)]
    walls=[label for label in allrows if label.startswith('taper_')]
    assert walls==[f'taper_{i}' for i in range(19,65)]
    archive=json.loads(source.read_bytes()); per=[];offsets={}
    for i,spec in enumerate(archive):
        f=Profile(spec);mid,slopes=f.root_clearance()
        def deficit(label,pose):
            u,h=allrows[label];return f.support(prepared['preimage'](u,*pose))-h
        other=min(deficit(label,pose).lower()for label,pose in zip(cap_rows,poses,strict=True)if label is not None)
        assert other>rational('0.0000009')
        qa={k:str(deficit(k,(240,-1)))for k in walls}
        for k,v in qa.items():
            val=arb(v).lower()
            if k not in offsets or val>arb(offsets[k]['lower']):offsets[k]={'lower':str(val),'profile':i}
        per.append(dict(index=i,parameters=list(map(str,f.x)),det=str(f.det),other_11_min=str(other),
                        root_clearance=mid,endpoint_derivative_clearance=slopes,qa=qa))
    return dict(complete=True,profile_count=len(per),input_sha256=sha256(source.read_bytes()).hexdigest(),
                host_sha256=sha256((comp/'ah-tapered-release.json').read_bytes()).hexdigest(),offsets=offsets,profiles=per)


if __name__=='__main__':
    a=ArgumentParser(description=__doc__)
    for name in ('computation','profiles','output'):a.add_argument(name,type=Path)
    args=a.parse_args()
    if not __debug__:raise SystemExit('Assertions must be enabled.')
    result=run(args.computation,args.profiles)
    result['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['profile_count'],sum(arb(t['lower'])>0 for t in result['offsets'].values()))
