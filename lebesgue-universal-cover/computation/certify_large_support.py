"""Certify an incremental area bound throughout a finite support cage."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from flint import arb, ctx

from certify_nine import ar
from certify_five_harmonic import old_common_polygon

SCALE=10**12
PI_UP=F(22,7)


def dot(u,v):
    return sum(a*b for a,b in zip(u,v,strict=True))


def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def normal(seed):
    s,q=seed
    q=F(q,10**8)
    assert s in (-1,1) and abs(q)<=1
    return s*(1-q*q)/(1+q*q),s*2*q/(1+q*q)


def outward(z):
    p=int(float(z)*SCALE)+2
    assert ar(F(p,SCALE))>z
    return F(p,SCALE)


def host_support_upper(params,u,pi,half=None):
    if half is None:
        half=[F(0)]*len(params)
    h=[F(1,2)]
    geometry=[(3,F(0),params[:2],F(0),half[:2])]
    geometry.extend((n,params[i],params[i+1:i+3],lip*half[i],half[i+1:i+3])
                    for (n,lip),i in zip(((5,F(163,1000)),(7,F(23,200)),(9,F(89,1000))),range(2,len(params),3)))
    for n,angle,t,angle_loss,translation_width in geometry:
        a=(pi/(2*n)).cos()
        r=1/(2*a)
        projections=[]
        for i in range(n):
            phi=ar(angle)+2*pi*i/n
            q=dot(tuple(map(ar,u)),(phi.cos(),phi.sin()))
            projections.append(r*q)
            if not q> -a:
                projections.append(1+r*q)
        h.append(max(outward(q) for q in projections)+dot(t,u)+angle_loss
                 +dot(translation_width,tuple(map(abs,u))))
    return max(h)


def new_support(u,t,terms):
    c=dot(tuple(map(ar,u)),(t.cos(),t.sin()))
    values=[arb(1),c]
    for _ in range(2,terms[-1][0]+1):
        values.append(2*c*values[-1]-values[-2])
    return ar(F(1,2))+sum(ar(coef)*values[k] for k,coef in terms)


def weights(us):
    a,b,c=us
    d=(cross(b,c),cross(c,a),cross(a,b))
    if max(d)<0:
        d=tuple(-x for x in d)
    assert min(d)>0
    d=tuple(x/sum(d) for x in d)
    assert tuple(sum(w*u[i] for w,u in zip(d,us,strict=True)) for i in (0,1))==(0,0)
    return d


def curvature(terms,pi):
    assert terms and all(k>=3 and k%2 for k,c in terms)
    assert list(terms)==sorted(terms) and len({k for k,c in terms})==len(terms)
    n=1024
    c=[(k,(1-k*k)*a) for k,a in terms]
    M=sum(k*k*abs(a) for k,a in c)
    target=F(1,100)
    pad=M*(2*PI_UP/n)**2/8
    for j in range(n):
        t=2*pi*j/n
        assert ar(F(1,2))+sum(ar(a)*(k*t).cos() for k,a in c)>ar(target+pad)
    return target


def generate(normals,base,terms,n):
    import numpy as np
    u=np.asarray(normals,dtype=float)
    i,j,k=np.array(list(combinations(range(len(u)),3))).T
    d=u[:,None,0]*u[None,:,1]-u[:,None,1]*u[None,:,0]
    W=np.array([d[j,k],d[k,i],d[i,j]])
    use=(W>1e-10).all(axis=0)|(W< -1e-10).all(axis=0)
    W=W[:,use]
    W/=W.sum(axis=0)
    i,j,k=i[use],j[use],k[use]
    ang=np.arctan2(u[:,1],u[:,0]); th=np.arange(n+1)*2*np.pi/n
    old=np.asarray(base,dtype=float)
    z=.5-old+sum(float(a)*np.cos(k*(ang[None,:]-th[:,None])) for k,a in terms)
    score = W[0]*z[:,i]+W[1]*z[:,j]+W[2]*z[:,k]
    picked=np.argmax(np.minimum(score[:-1],score[1:]),axis=1)
    return np.c_[i[picked],j[picked],k[picked]].tolist()


def verify(normals,base,terms,grid,rows,pi):
    assert len(rows)==grid
    M=sum(k*k*abs(a) for k,a in terms)
    remainder=M*(2*PI_UP/grid)**2/8
    lower=None
    used=set()
    previous={}
    for j,indices in enumerate(rows):
        assert len(set(indices))==3 and all(type(i) is int and 0<=i<len(normals) for i in indices)
        used.update(indices)
        us=[normals[i] for i in indices]
        w=weights(us)
        now={i:new_support(normals[i],2*pi*(j+1)/grid,terms) for i in indices}
        for i in indices:
            if i not in previous:
                previous[i]=new_support(normals[i],2*pi*j/grid,terms)
        for sample in (previous,now):
            e=sum(ar(a)*(sample[i]-ar(base[i])) for a,i in zip(w,indices,strict=True))
            bound=-outward(-e)-remainder
            lower=bound if lower is None else min(lower,bound)
        previous=now
    return lower,sorted(used)


def containment(spec,terms,pi):
    assert [[k,str(a)] for k,a in terms]==spec['terms']
    params=list(map(F,spec['params']))
    rho,tx,ty=map(F,spec['test_pose'])
    hull,_=old_common_polygon(params,[F(0)]*len(params),pi,384)
    worst=None
    for p,q in zip(hull,hull[1:]+hull[:1],strict=True):
        x,y=q[1]-p[1],p[0]-q[0]
        length=arb(x*x+y*y).sqrt()
        angle=arb.atan2(arb(y),arb(x))-ar(rho)
        h=ar(F(1,2))+sum(ar(a)*(k*angle).cos() for k,a in terms)
        slack=arb(x*p[0]+y*p[1])/(length*SCALE)-h-ar(x*tx+y*ty)/length
        s=-outward(-slack)
        assert s>0
        worst=s if worst is None else min(worst,s)
    return str(worst)


def cap_required(h, area):
    # Cap depth in the first grid interval attaining the given area.
    d=10**10
    lo,hi=0,d
    while hi-lo>1:
        mid=(lo+hi)//2
        q=F(mid,d)
        if q**4/4>area**2*((h+q)**2-F(1,4)):
            hi=mid
        else:
            lo=mid
    q=F(hi,d)
    assert q**4/4>area**2*((h+q)**2-F(1,4))
    return q


def stationary_drift(center, balls):
    change=[outward(abs(arb(v)-ar(c))) for v,c in zip(balls,center,strict=True)]
    delta=[ar(change[0]**2+change[1]**2).sqrt()]
    for i,l in zip(range(2,len(center),3),(F(163,1000),F(23,200),F(89,1000))):
        delta.append(ar(l*change[i])+ar(change[i+1]**2+change[i+2]**2).sqrt())
    return max(outward(d) for d in delta)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--generate',type=Path)
    g.add_argument('--verify',type=Path)
    p.add_argument('--containment',type=Path)
    a=p.parse_args()
    ctx.prec=160
    pi=arb.pi()
    assert pi<ar(PI_UP)
    source=a.input.read_bytes()
    spec=json.loads(source)
    coeffs=[(k,F(v)) for k,v in spec['terms']]
    rmin=curvature(coeffs,pi)
    centers=[list(map(F,v)) for v in spec['centers']]
    assert [len(c) for c in centers]==[8,11]
    drifts=[stationary_drift(c,v['placement']) for c,v in zip(centers,spec['stationary_candidates'],strict=True)]
    normals=[normal(seed) for seed in spec['normal_seeds']]
    assert all(dot(u,u)==1 for u in normals)
    base=[max(host_support_upper(params,u,pi) for params in centers) for u in normals]
    stamp=sha256(source).hexdigest()
    grid=spec['angle_cells']
    assert all(1/(2*(pi/(2*j)).cos())<ar(F(3,5)) for j in (3,5,7,9))
    assert all((pi/(2*j)).tan()/2<ar(l) for j,l in ((5,F(163,1000)),(7,F(23,200)),(9,F(89,1000))))
    err=F(3,5)*max((centers[0][j]-centers[1][j])**2 for j in (2,5))/8
    if a.generate:
        table=[]
    else:
        w=json.loads(a.verify.read_bytes())
        assert w['input_sha256']==stamp
        table=w['indices']
        assert len(table)==len(spec['area_guards'])+len(spec.get('traced_boxes',[]))
    guards=[]
    for j,row in enumerate(spec['area_guards']):
        padding,halfwidth,area=map(F,row)
        assert padding>0 and halfwidth>0 and area>0
        assert max(drifts)<padding
        assert ar(halfwidth)*(arb(2).sqrt()+ar(F(163,1000)))+ar(err)<ar(padding)
        supports=[h+padding for h in base]
        expanded=[b+cap_required(b,area) for b in supports]
        if a.generate:
            table.append(generate(normals,expanded,coeffs,grid))
        slack,used=verify(normals,expanded,coeffs,grid,table[j],pi)
        assert slack>0
        guards.append({'cage_padding':str(padding),'raw_coordinate_halfwidth':str(halfwidth),
                       'strict_area_increase':str(area),'balanced_cap_slack':str(slack),
                       'normals_used':len(used),'cage':[[spec['normal_seeds'][i],str(supports[i])] for i in used]})
    out={'input_sha256':stamp,'curvature_lower':str(rmin),'angular_cells':grid,
         'tube_interpolation_padding':str(err),'stationary_support_drifts':list(map(str,drifts)),
         'area_guards':guards}
    trace=[]
    for j,s in enumerate(spec.get('traced_boxes',[]),start=len(guards)):
        lo,hi=[list(map(F,s[k])) for k in ('lo','hi')]
        center=[(a+b)/2 for a,b in zip(lo,hi,strict=True)]
        half=[(b-a)/2 for a,b in zip(lo,hi,strict=True)]
        assert len(center)==11 and min(half)>0
        supports=[host_support_upper(center,u,pi,half) for u in normals]
        area=F(s['strict_area_increase'])
        assert area>0
        expanded=[b+cap_required(b,area) for b in supports]
        if a.generate:
            table.append(generate(normals,expanded,coeffs,grid))
        slack,used=verify(normals,expanded,coeffs,grid,table[j],pi)
        assert slack>0, (s['part'],s['pos'],float(slack))
        trace.append({'part':s['part'],'pos':s['pos'],'normals_used':len(used),
                      'strict_area_increase':str(area),'balanced_cap_slack':str(slack)})
    if trace:
        out['traced_box_area_guards']=trace
    if a.containment:
        out['explicit_zero_gain_configuration_padding']=containment(json.loads(a.containment.read_bytes()),coeffs,pi)
    if a.generate:
        a.generate.write_text(json.dumps({'input_sha256':stamp,'indices':table},separators=(',',':'))+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
