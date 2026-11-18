"""Construct a marked Reuleaux pentagon in the Gibbs hexagon and evaluate its genuine support."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from flint import ctx
from certify_sigma_extended_graph import geometry


@lru_cache(maxsize=2)
def environment(archive,sigma):
    precision=ctx.prec
    try:
        ns,_,_=geometry(Path(archive).resolve(),[(F('.429'),F('.009'))],sigma=str(sigma),far=F('.273'),negative=F('.306'),eps=F('.0102'))
        return ns
    finally:ctx.prec=precision


def body(ns,w,deformation):
    A,add,sub,mul,dot,reflect=(ns[k] for k in ('A','add','sub','mul','dot','reflection'));l,t=ns['chart'](add(ns['z'],mul(A(w),ns['D'])));p=ns['psi'](l);t-=A(deformation)
    vertices=[add(ns['l0'],mul(l,ns['el'])),reflect(add(ns['l0'],mul(p,ns['el']))),reflect(ns['X'](p,-t)),add(ns['mid'],mul(t,ns['en'])),ns['X'](l,t)]
    edges=[sub(vertices[i],vertices[j]) for i,j in ns['orders']];assert all(dot(e,e).contains(1) for e in edges);assert all(ns['cross'](edges[i],edges[(i+1)%10])>0 for i in range(10))
    assert all(dot(d:=sub(vertices[i],vertices[j]),d)<1 for i in range(5) for j in range(i) if (i,j) not in ns['orders'])
    def h(u):
        if u[0].is_zero() and u[1].is_zero():return A(0)
        radius=ns['norm'](u);direction=mul(1/radius,u)
        return radius*ns['supports'](vertices,edges,(A(1),A(0)),direction)
    return vertices,h


def support(archive,sigma,w,deformation):
    return body(environment(str(Path(archive).resolve()),str(sigma)),F(w),F(deformation))[1]
