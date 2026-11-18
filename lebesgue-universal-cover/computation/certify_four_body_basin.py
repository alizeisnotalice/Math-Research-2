"""Explicit four-body basin via correlated Taylor models and seven frozen fans."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import flint
from flint import arb,arb_mat,ctx
from four_contact_taylor import Model
from certify_four_body_contact_stationary import (
    BOX,PIECE_DATA,TIES,add,chain,det,dot,geometry,interval,left,normal,pi,primitive_key,q,r,regular,scale,sub,unit_root,
)
from certify_four_body_contact_transverse import EVENT_ORDER,SIMPLEX,connector,face_orders,plus,rigid_velocities

CENTERS=('0.51250444532618500354392','1.577473015380594982391575')
RAW_CENTER=tuple(map(F,('-0.00255163398964655','-0.00193331458560072',
    '0.02867263048929559','0.02508688430651930','-0.00328479589201898',
    '0.70156070130903857','0.00288748433564455','-0.01271885155848397')))
RAW_WIDTH=F('0.000001')
H_MIN=F('0.1')
CUBE=F('0.0001')
CHOICES=list(product(range(4),range(4),range(4),range(4),range(2),range(2)))


def rotate(p,co,si):
    return add(scale(co,p),scale(si,left(p)))


def area_model(pieces):
    start=Model();total=Model()
    for (_,stop,_,rad,C),(_,_,_,rad2,C2) in zip(pieces,pieces[1:]+pieces[:1]):
        n0,n1=normal(start),normal(stop)
        total+=rad*det(C,sub(n1,n0))+rad**2*(stop-start)+det(add(C,scale(rad,n1)),add(C2,scale(rad2,n1)))
        start=stop
    return total/2


def mapped_data(data,x,y):
    def ev(model):
        return model.evaluate(x,y) if isinstance(model,Model) else model
    def pt(p):
        return tuple(ev(v) for v in p)
    pieces,t,ns,p=data
    return [(tag,ev(stop),i,rad,pt(c)) for tag,stop,i,rad,c in pieces],pt(t),{k:pt(v) for k,v in ns.items()},{n:[pt(v) for v in Vs] for n,Vs in p.items()}


def topology(data):
    pieces,t,ns,p=data
    assert all(ns[k][0]>0 for k in ('C','D','E','F'))
    assert all(ns[k][0]<0 for k in ('X','Y','Z'))
    cuts={'start':arb(0),'end':2*pi}
    rotations={3:arb(0),5:t[2],7:t[5]}
    for n,f in rotations.items():
        for i in range(n):
            for end in (-1,1):
                angle=f+2*pi*i/n+end*pi/(2*n)
                if angle<0: angle+=2*pi
                if angle>2*pi: angle-=2*pi
                assert 0<angle<2*pi
                cuts[n,i,end]=angle
    for k,angle,*_ in pieces:
        if k not in cuts:cuts[k]=angle
    keys=sorted(cuts,key=lambda k:float(cuts[k].mid()))
    assert (keys[0],keys[-1])==('start','end') and all(cuts[k]<cuts[j] for k,j in zip(keys,keys[1:]))
    active=0
    for k,j in zip(keys,keys[1:]):
        l,h=cuts[k],cuts[j];assert h-l<pi
        owner,entry,rad=PIECE_DATA[active]
        mid=(l+h)/2;nl,nh=normal(l),normal(h)
        if owner:assert primitive_key(owner,mid,rotations[owner])==(owner,entry,rad)
        c=p[owner][entry] if owner else (0,0)
        for n in (0,3,5,7):
            if n==owner:continue
            if n:
                _,i,rr=primitive_key(n,mid,rotations[n]);d=p[n][i]
            else:rr,d=r,(0,0)
            v=sub(c,d)
            for tag,endpoint in ((k,nl),(j,nh)):
                if {n,owner}<=TIES.get(tag,set()):continue
                assert rad-rr+dot(v,endpoint)>0,('support_endpoint',tag,owner,n,rad-rr+dot(v,endpoint))
            nv=scale(-1,v)
            if not (det(nl,nv)<0 or det(nv,nh)<0):
                assert rad-rr-dot(v,v).sqrt()>0,('support_interior',k,j,owner,n)
        if j==pieces[active][0]:active+=1
    assert active==len(pieces)


def fan_parts(data,correlated=False):
    pieces,t,ns,p=data
    if correlated:
        a,b,c,d=(ns[k] for k in ('A','B','C','D'))
        orders={'A':[(3,p[3][0]),(5,add(p[5][3],a)),(0,scale(r,a)),(7,p[7][0])],
                'B':[(5,p[5][1]),(0,scale(r,b)),(7,p[7][1]),(3,p[3][1])],
                'C':[(5,add(p[5][2],c)),(7,add(p[7][2],c)),(0,scale(r,c))],
                'D':[(0,scale(r,d)),(7,add(p[7][2],d)),(3,p[3][0])]}
        for tag in ('A','B'):
            orders[tag+'opp']=[(i,sub(pt,ns[tag])) for i,pt in reversed(orders[tag])]
    else:
        orders=face_orders(p,ns)
        assert [[i for i,_ in orders[k]] for k in EVENT_ORDER]==[[3,5,0,7],[7,0,5,3],[5,0,7,3],[3,7,0,5],[5,7,0],[0,7,3]]
    actual=[];arcs=[]
    start=arb(0);arc_g=[arb(0)]*6
    for (tag,stop,i,rad,C),(_,_,j,rad2,C2) in zip(pieces,pieces[1:]+pieces[:1]):
        n0,n1=normal(start),normal(stop)
        P=add(C,scale(rad,n1));Q=add(C2,scale(rad2,n1))
        if i in (5,7) and rad==1:
            diff=sub(n1,n0);arcs.append((i,rad,diff))
            arc_g=plus(arc_g,[rad*det(v,diff) for v in rigid_velocities(i,C,t)])
        if i!=j:
            if tag in orders:
                seq=orders[tag];assert seq[0][0]==i and seq[-1][0]==j
                actual.append((tag,seq))
            else:
                actual.append((tag,[(i,P),(j,Q)]))
        start=stop
    for index in SIMPLEX:
        opts=dict(zip(EVENT_ORDER,CHOICES[index]));links=[];g=list(arc_g)
        for tag,seq in actual:
            if tag in opts:
                k=opts[tag];mid=len(seq)-2
                keep=[seq[0],*(item for j,item in enumerate(seq[1:-1]) if k&(1<<(mid-j-1))),seq[-1]]
            else:keep=seq
            for old,next_ in zip(keep,keep[1:]):
                links.append((old,next_));g=plus(g,connector(*old,*next_,t))
        yield g,links,arcs,t


def bound_fan_hessian(links,arcs,t,D):
    """Six rigid-motion variables, all with absolute value at most D."""
    delta=arb(0,D); co,si=delta.cos(),delta.sin();dh=[[(0,0)]*6 for _ in range(6)]
    trans={5:add(t[3:5],(delta,delta)),7:add(t[6:8],(delta,delta))}
    zero=(0,0); unit5=(1,0);unit7=(0,1)
    H=[[arb(0) for _ in range(6)] for _ in range(6)]
    for i,rad,diff in arcs:
        j=0 if i==5 else 3
        rotated=rotate(diff,co,si)
        H[j][j]-=rad*det(trans[i],rotated)/2
        for k,e in enumerate((unit5,unit7),j+1):
            H[j][k]+=rad*det(e,left(rotated))/2;H[k][j]=H[j][k]
    def moved(i,at):
        first={};second={}
        if i in (0,3):return at,first,second
        j=0 if i==5 else 3
        v=rotate(sub(at,t[j+3:j+5]),co,si)
        first={j:left(v),j+1:unit5,j+2:unit7};second={j:scale(-1,v)}
        return add(trans[i],v),first,second
    min_orientation=None
    for (i,P),(j,Q) in links:
        V,dv,hv=moved(i,P);W,dw,hw=moved(j,Q)
        orientation=det(V,W).lower()
        min_orientation=orientation if min_orientation is None else min(min_orientation,orientation)
        for c,dp in hv.items():H[c][c]+=det(dp,W)/2
        for c,dq in hw.items():H[c][c]+=det(V,dq)/2
        for c,dp in dv.items():
            for d,dq in dw.items():
                increment=det(dp,dq)/2
                H[c][d]+=increment;H[d][c]+=increment
    norm=max(sum(abs(a).abs_upper() for a in row) for row in H)
    return min_orientation,norm



def prepare_gradient(data):
    fans=list(fan_parts(data,correlated=True))
    G=[[Model(1)]*7,*[[f[0][j] for f in fans] for j in range(6)]]
    z=arb(0)
    inverse=arb_mat([[entry.evaluate(z,z) for entry in row] for row in G]).inv()
    B=[[arb(round(float(inverse[i,j].mid())*10**12))/10**12 for j in range(7)] for i in range(7)]
    Bmat=arb_mat(B)
    det_B=Bmat.det()
    assert det_B>0 or det_B<0
    P=[[sum(B[i][k]*G[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
    return P,Bmat


def simplex_margin(polynomial_matrix,x,y):
    P,B=polynomial_matrix
    matrix=arb_mat([[v.evaluate(x,y) for v in row] for row in P])
    inverse=matrix.inv()*B
    cube=q(CUBE)
    return min((inverse[i,0]-cube*sum(abs(inverse[i,j]) for j in range(1,7))).lower() for i in range(7))


def cell(data,value,correlated,R,x,y,trans_radius):
    H=[[value.second(i,j,x,y,R) for j in range(2)] for i in range(2)]
    H[0][0]-=q(H_MIN);H[1][1]-=q(H_MIN)
    minor=H[0][0].lower()
    determinant=det(*H).lower()
    assert minor>0 and determinant>0,('tangent_curvature',minor,determinant)
    local=mapped_data(data,x,y)
    topology(local)
    fans=list(fan_parts(local))
    cube_margin=simplex_margin(correlated,x,y)
    assert cube_margin>0,('transverse_cube',cube_margin)
    h=[bound_fan_hessian(*f[1:],trans_radius) for f in fans]
    orientation=min(z[0] for z in h); remainder=max(z[1] for z in h)*trans_radius/2
    linear=(q(CUBE)-remainder).lower()
    assert orientation>0 and linear>0,('fan_remainder',orientation,linear)
    return {'tangent_minor':minor,'tangent_determinant':determinant,'gradient_cube_margin':cube_margin,
            'connector_orientation':orientation,'transverse_linear_coefficient':linear}


def mapped_raw(data,R,A,B,tr,raw_width):
    raw=[arb(q(v),q(raw_width)) for v in RAW_CENTER]
    t3=raw[:2]
    aa=unit_root(add(t3,regular(3,0)),r,1); bb=unit_root(add(t3,regular(3,1)),r,-1)
    assert aa[0]>0 and bb[1]>0
    ux=((aa[1]/aa[0]).atan()-A)/R
    uy=(pi/2-(bb[0]/bb[1]).atan()-B)/R
    assert abs(ux)<q('.5') and abs(uy)<q('.5'),(ux,uy)
    rest=mapped_data(data,ux,uy)[1][2:]
    result=[abs(a-b).abs_upper() for a,b in zip(raw[2:],rest)]
    assert max(result)<tr,result
    return [str(ux),str(uy)],list(map(str,result))


def main():
    parser=ArgumentParser()
    parser.add_argument('--radius',type=F,default=F('.00004'))
    parser.add_argument('--transverse',type=F,default=F('.0001'))
    parser.add_argument('--degree',type=int,default=12)
    parser.add_argument('--depth',type=int,default=16)
    parser.add_argument('--raw-width',type=F,default=RAW_WIDTH)
    args=parser.parse_args()
    Model.degree=args.degree
    A,B=map(q,CENTERS);R=q(args.radius)
    alpha=Model({(0,0):A,(1,0):R});beta=Model({(0,0):B,(0,1):R})
    data=chain(alpha,beta)
    V=area_model(data[0]);correlated=prepare_gradient(data);tr=q(args.transverse)
    assert all(abs((interval(*bounds)-c)/R)<q('.5') for bounds,c in zip(BOX,(A,B)))
    raw_angles,raw_trans=mapped_raw(data,R,A,B,tr,args.raw_width)
    floors={};counts=Counter(); deepest=0
    def visit(lo,hi,depth):
        nonlocal deepest
        counts['visited']+=1;deepest=max(deepest,depth)
        x,y=(interval(a,b) for a,b in zip(lo,hi))
        try:
            values=cell(data,V,correlated,R,x,y,tr)
        except (AssertionError,ValueError,ZeroDivisionError) as e:
            if depth==args.depth:raise
            counts[str(e.args[0][0]) if e.args and isinstance(e.args[0],tuple) else 'other']+=1
            j=depth%2;m=(lo[j]+hi[j])/2
            visit(lo,hi[:j]+(m,)+hi[j+1:],depth+1)
            visit(lo[:j]+(m,)+lo[j+1:],hi,depth+1);return
        counts['terminal']+=1
        for tag,val in values.items():floors[tag]=min(floors.get(tag,val),val)
    half=F(1,2);visit((-half,-half),(half,half),0)
    assert counts['visited']==2*counts['terminal']-1
    source=Path(__file__)
    print(json.dumps({'python':platform.python_version(),'flint':flint.__version__,'precision':ctx.prec,
       'script_sha256':sha256(source.read_bytes()).hexdigest(),
       'taylor_script_sha256':sha256(source.with_name('four_contact_taylor.py').read_bytes()).hexdigest(),
       'taylor_degree':args.degree,'taylor_polydisk_radius':str(args.radius),'area_taylor_remainder_upper':str(V.error),
       'tangent_center':CENTERS,'tangent_halfwidth':str(args.radius/2),'tangent_curvature_coefficient_lower':str(H_MIN),
       'transverse_halfwidth':str(args.transverse),'gradient_cube_radius':str(CUBE),'counts':dict(counts),
       'deepest':deepest,'strict_margins':{k:str(v) for k,v in floors.items()},
       'raw_center':list(map(str,RAW_CENTER)),'raw_halfwidth':str(args.raw_width),'raw_normalized_tangent':raw_angles,
       'raw_transverse_sup':raw_trans},indent=2))


if __name__=='__main__':
    main()
