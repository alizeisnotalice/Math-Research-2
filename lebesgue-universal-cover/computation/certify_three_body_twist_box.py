"""Certified coercivity on a near-zero angular box for D, B3 and B5."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, arb_series, ctx
from certify_three_body_symmetric_slice import R3, add, c36, c72, det, dot, g, interval, pi, q, r, scale, sub
from certify_three_body_twist_jet import J, rotate, other_root, argument_right

ctx.cap = 3


def contacts(a, y, phi, upper):
    right=(R3+a,y)
    ux=a-R3/2
    up3,dn3=(ux,r+y),(ux,-r+y)
    length=dot(right,right).sqrt()
    K=(dot(right,right)-r**2).sqrt()
    e=scale(1/length,right)
    f=J(e)
    f5=(-phi.sin(),phi.cos())
    d,o=dot(f5,f),dot(f5,e)
    plus=add(scale(r/length,e),scale(K/length,f))
    minus=add(scale(r/length,e),scale(-K/length,f))
    w0=add(scale(-length+2*g*K*d,e),scale(g*r*o/K,f))
    t5=add(w0,rotate((c36,0),phi))
    w1,w2=sub(w0,scale(g,f5)),add(w0,scale(g,f5))
    up5=add(w0,rotate((c36+c72,r),phi))
    dn5=add(w0,rotate((c36+c72,-r),phi))
    m=1 if upper else -1
    t=y/(ux-(ux**2+(y+m)*y).sqrt())
    top_n=(-2*t/(1+t**2),(1-t**2)/(1+t**2))
    return dot(top_n,up5)-r, (right,up3,dn3,up5,dn5,w1,w2,plus,minus,top_n,t5)


def chain(data, phi, upper):
    right,up3,dn3,up5,dn5,w1,w2,plus,minus,top_n,_ = data
    gp,gm = argument_right(plus),argument_right(minus)
    dp,dm = argument_right(other_root(w1,plus)),argument_right(other_root(w2,minus))
    top = pi/2-(top_n[0]/top_n[1]).atan()
    zero=0*phi
    origin=(zero,zero)
    parts = [(gp,0,right,3), (dp,r,origin,0), (3*pi/10+phi,1,w1,5), (top,0,up5,5)]
    if not upper:
        parts.append((pi/2,1,dn3,3))
    parts.extend([(5*pi/6,0,up3,3), (pi+gm,1,right,3), (pi+gp,r,origin,0), (7*pi/6,1,right,3)])
    if upper:
        parts.extend([(3*pi/2,0,dn3,3), (pi+top,1,up3,3)])
    else:
        parts.append((pi+top,0,dn3,3))
    parts.extend([(3*pi/2+phi,1,up5,5), (17*pi/10+phi,0,dn5,5), (2*pi+dm,1,w2,5),
                  (2*pi+gm,r,origin,0), (2*pi,0,right,3)])
    return parts,(gp,dp,gm,dm,top)


def parts_eval(parts,zero):
    first=zero
    for i,(last,rad,center,owner) in enumerate(parts):
        _,rad2,center2,owner2=parts[(i+1)%len(parts)]
        n1=(first.cos(),first.sin())
        n2=(last.cos(),last.sin())
        yield (first,last,rad,center,owner,n1,n2,add(center,scale(rad,n2)),add(center2,scale(rad2,n2)),owner2)
        first=last


def area(parts,zero):
    total=0*zero
    for x,y,rad,ctr,_,n1,n2,P,Q,_ in parts_eval(parts,zero):
        total += rad*det(ctr,sub(n2,n1)) + rad**2*(y-x) + det(P,Q)
    return total/2


def short_area(details,phi,upper):
    p,u,v,top5,_,wplus,wminus,_,_,nT,_=details
    length=dot(p,p).sqrt()
    K=(dot(p,p)-r**2).sqrt()
    rotated_parallel=p[0]*phi.cos()+p[1]*phi.sin()
    d=rotated_parallel/length
    o=(-p[0]*phi.sin()+p[1]*phi.cos())/length
    T=pi/2-(nT[0]/nT[1]).atan()
    A=14*pi/15-2*r**2*(K/r).atan()+2*(1-r**2)*((wminus[1]/wminus[0]).atan()-(wplus[1]/wplus[0]).atan())+phi
    if not upper:
        A += pi-2*T
    U=u if upper else v
    D=K+4*g*rotated_parallel+4*g**2*K*d**2-2*g**2*r*o**2/K+2*g*(c36+c72)-(pi/6).cos()-(3*pi/10).cos()
    return (A+D+2*det(sub(top5,U),nT))/2,T


def jets(a,y,phi,upper):
    s=arb_series([phi,1]); t=arb_series([y,1])
    E_y=contacts(a,t,arb_series([phi]),upper)[0][1]
    E_phi=contacts(a,arb_series([y]),s,upper)[0][1]
    v=-E_phi/E_y
    e2=contacts(a,arb_series([y,v]),s,upper)[0][2]
    w=-e2/E_y
    details = contacts(a,arb_series([y,v,w]),s,upper)[1]
    H,T=short_area(details,s,upper)
    return E_y, v, H[2], T[1]


def solve_interval(a,phi,upper, full_y):
    """Use a strict point range on E_y and a midpoint residual to enclose roots."""
    ey=contacts(a,arb_series([full_y,1]),arb_series([phi]),upper)[0][1]
    assert ey > q('.7'), ey
    ac,pc = a.mid(),phi.mid()
    guess = arb(0)
    for _ in range(3):
        effect=contacts(ac,arb_series([guess,1]),arb_series([pc]),upper)[0]
        guess = (guess-effect[0]/effect[1]).mid()
    assert full_y.contains(guess)
    residual0=contacts(ac,guess,pc,upper)[0]
    grad_a=contacts(arb_series([a,1]),arb_series([guess]),arb_series([phi]),upper)[0][1]
    grad_phi=contacts(arb_series([a]),arb_series([guess]),arb_series([phi,1]),upper)[0][1]
    residual=residual0+grad_a*(a-ac)+grad_phi*(phi-pc)
    assert abs(residual) < ey.lower()*(full_y.upper()-abs(guess)), (residual,ey,guess)
    root=guess-residual/ey
    assert full_y.contains(root), (full_y,root)
    return root,ey


# Transverse variables in order y3, b, z.
BASIS = [{3:(0,1),5:(0,0)}, {3:(0,0),5:(1,0)}, {3:(0,0),5:(0,1)}]


def velocity(basis,owner):
    return basis.get(owner,(0,0))


def gradient(parts):
    values = [arb(0),arb(0),arb(0)]
    for _,_,rad,_,owner,n1,n2,P,Q,owner2 in parts_eval(parts,arb(0)):
        for j, basis in enumerate(BASIS):
            V,W = velocity(basis,owner),velocity(basis,owner2)
            values[j] += (rad*det(V,sub(n2,n1)) + det(V,Q) + det(P,W))/2
    return values


def triangle_gradient(P,M,Q,i,j,k):
    values=[]
    for basis in BASIS:
        V,W,Z = (velocity(basis,owner) for owner in (i,j,k))
        values.append((det(sub(W,V),sub(Q,P)) + det(sub(M,P),sub(Z,V)))/2)
    return values


def det3(m):
    return (m[0][0]*det(m[1][1:],m[2][1:]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
             +m[0][2]*det(m[1][:2],m[2][:2]))


def transpose(m):
    return list(zip(*m))


def solve3(m,b):
    denominator=det3(m)
    columns=transpose(m)
    result=[]
    for j in range(3):
        replaced=[*columns]
        replaced[j]=b
        result.append(det3(transpose(replaced))/denominator)
    return result


def transverse_bounds(details,parts,upper):
    P,u,v,top5,_,w1,w2,np,nm,nT,_=details
    top3=u if upper else add(v,nT)
    T1=triangle_gradient(P,add(w1,np),scale(r,np),3,5,0)
    T2=triangle_gradient(scale(r,nm),add(w2,nm),P,0,5,3)
    T3=triangle_gradient(top5,scale(r,nT),top3,5,0,3)
    matrix=[T1,T2,T3]
    ell=gradient(parts)
    weights=solve3(transpose(matrix),ell)
    inverses=transpose([solve3(matrix,[int(i==j) for i in range(3)]) for j in range(3)])
    inverse_max=max(abs(x).abs_upper() for row in inverses for x in row)
    return weights, 1/inverse_max, matrix



def base_topology(a):
    b=contacts(a,arb(0),arb(0),True)[1][-1][0]
    right_B3,left_B3=R3+a,a-R3/2
    from certify_three_body_symmetric_slice import R5
    right_B5=R5+b
    w=(b-c36,-g)
    e18,e30,e54=[((pi*k/180).cos(),(pi*k/180).sin()) for k in (18,30,54)]
    return min(t.lower() for t in [
        right_B3-right_B5,
        right_B5*e18[0]-r,
        r-right_B3*e30[0],
        -r-dot(w,e30),
        1-dot((right_B3-w[0],g),e30),
        w[0]*e30[1]+g*e30[0],
        (right_B3-w[0])*e18[1]-g*e18[0],
        -w[0],g,left_B3-w[0],r-g,
        -dot((left_B3-w[0],-r+g),e30),
        r+dot(w,e54),
        dot((c72+b,r),e54)-r,
    ])


def shifted_bound(details,parts,actual,upper):
    P,u,v,top5,_,w1,w2,np,nm,nT,t5=details
    top3=u if upper else add(v,nT)
    shifts={3:(arb(0),actual[0]-P[1]),5:(actual[1]-t5[0],actual[2]-t5[1]),0:(0,0)}
    shifts={j:tuple(arb(0,abs(value).abs_upper()) if isinstance(value,arb) else arb(0) for value in row) for j,row in shifts.items()}
    cl=[]
    def check(A,B,i,j):
        cl.append(det(add(A,shifts[i]),add(B,shifts[j])))
    for _,_,_,_,i,_,_,C,D,j in parts_eval(parts,arb(0)):
        if i!=j:
            check(C,D,i,j)
    for n,C,D,E,i,j,k in [
        (np,P,add(w1,np),scale(r,np),3,5,0),
        (nm,scale(r,nm),add(w2,nm),P,0,5,3),
        (nT,top5,scale(r,nT),top3,5,0,3),
    ]:
        check(C,D,i,j)
        check(D,E,j,k)
        check(sub(E,n),sub(D,n),k,j)
        check(sub(D,n),sub(C,n),j,i)
    return min(x.lower() for x in cl), min(abs(shifts[j][u]).abs_upper() for j,u in [(3,1),(5,0)])


def certified_cell(a,phi,upper,actual,full_y):
    y,ey0=solve_interval(a,phi,upper,full_y)
    ey,_,coef,T1=jets(a,y,phi,upper)
    _,detail=contacts(a,y,phi,upper)
    _,_,_,_,_,wplus,wminus,nplus,nminus,nT,_=detail
    branch_margin=min(x.lower() for x in [nplus[0],nminus[0],other_root(wplus,nplus)[0],
        other_root(wminus,nminus)[0],nT[1],-wplus[0],-wminus[0]])
    parts,angles=chain(detail,phi,upper)
    gp,dp,gm,dm,T=angles
    layout=min(x.lower() for x in [gp-pi/10-phi,pi/6-gp,gm+pi/6,-pi/10+phi-gm,
        dp-pi/6,3*pi/10+phi-dp,dm+3*pi/10-phi,-pi/6-dm,T-pi*q('.49'),pi*q('.51')-T])
    weights,c0,_=transverse_bounds(detail,parts,upper)
    w=max(abs(v).abs_upper() for v in weights)
    lc,motion=shifted_bound(detail,parts,actual,upper)
    error=abs(weights[1]-weights[0]-weights[2]).abs_upper()/2
    inverse_gap=c0-3*motion/2
    sharp=(1-w)*inverse_gap-error*motion
    return {'root_derivative':ey0.lower(),'angular_curvature_half':coef.lower(),
            'upper_top_slope_margin':(1-abs(T1)).lower(), 'layout':layout,
            'argument_branch_margin':branch_margin, 'fan_connector_determinant':lc,
            'weight_margin':(1-w).lower(), 'inverse_gap':inverse_gap.lower(), 'sharp_transverse_coefficient':sharp.lower()}


TARGETS = {
    'hardest': (('-.01362','-.01323'), '0.002455', (('-.00038','.00001'),('.02230','.02269'),('-.00001','.000385'))),
    'near': (('-.0145','-.0094'), '0.0045', (('-.00076','.00076'),('.0219','.02346'),('-.00154','.000385'))),
    'wide': (('-.017','-.0076'), '0.0025', (('-.0045','.0045'),('.0185','.0274'),('-.0045','.0045'))),
}


def subintervals(first,last,count):
    first,last=F(first),F(last)
    for j in range(count):
        yield interval(first+(last-first)*F(j,count),first+(last-first)*F(j+1,count))


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('target',choices=TARGETS)
    parser.add_argument('--a-count',type=int,required=True)
    parser.add_argument('--phi-count',type=int,required=True)
    args=parser.parse_args()
    aa,phi_end,bbox=TARGETS[args.target]
    actual=tuple(interval(*pair) for pair in bbox)
    full_y=interval('-.001','.001')
    seen={}
    worst={}
    phi_items=list(subintervals(0,phi_end,args.phi_count))
    for i,a in enumerate(subintervals(*aa,args.a_count)):
        b0=base_topology(a)
        assert b0>0,('zero_angle_layout',i,b0)
        for j,phi in enumerate(phi_items):
            for upper in (False,True):
                values=certified_cell(a,phi,upper,actual,full_y)
                for name,val in values.items():
                    assert val>0,(name,i,j,upper,val,a,phi)
                    if name not in seen or val<seen[name]:
                        seen[name]=val;worst[name]=(i,j,upper)
    result={'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'target':args.target,
        'box_a':aa,'box_phi':['0',phi_end],'box_y_b_z':bbox,
        'partition':[args.a_count,args.phi_count],'both_analytic_branches':True,
        'evaluations':2*args.a_count*args.phi_count,'strict_margins':{k:str(v) for k,v in seen.items()},
        'worst_index':worst}
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
