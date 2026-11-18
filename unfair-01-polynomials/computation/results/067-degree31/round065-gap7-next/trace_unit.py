"""Finite-field Buchberger unit witness with deferred multiplier expansion.

Each stored reduction has edges to earlier polynomials. The terminal 1 is
expanded only once. The returned identity is always multiplied back.
"""
from sympy.polys.orderings import grevlex


def unit_witness(initial):
    if not initial:
        return None
    original=initial[0].ring
    R=original.clone(order=grevlex)
    polys=[p.set_ring(R) for p in initial]
    k=len(polys)
    edges=[None]*k
    basis=set()
    pairs=set()
    origin={}
    order=R.order
    lcm=R.monomial_lcm
    quotient=R.monomial_div

    def lift(index):
        weights={index:R.one}
        for i in range(index,k-1,-1):
            v=weights.pop(i,R.zero)
            if not v:continue
            for j,u in edges[i]:
                weights[j]=weights.get(j,R.zero)+v*u
        result=[weights.get(i,R.zero).set_ring(original) for i in range(k)]
        assert sum(a*b for a,b in zip(initial,result,strict=True))==original.one
        return result

    def normalize(f,source):
        if not f:return None
        reducers=sorted(basis,key=lambda j:order(polys[j].LM))
        if reducers:
            q,r=f.div([polys[j] for j in reducers])
        else:
            q,r=[],f
        if not r:return None
        inv=R.domain.one/r.LC
        r=r*inv
        if r in origin:
            return None
        # Each edge is a sparse factor times an earlier node.
        incoming=[(i,g*inv) for i,g in source if g]
        incoming.extend((j,-g*inv) for j,g in zip(reducers,q,strict=True) if g)
        i=len(polys)
        origin[r]=i
        polys.append(r)
        edges.append(incoming)
        return i

    def install(i):
        # Apply the relatively prime leading-monomial criterion and the
        # leading-lcm criterion to new pairs, then remove obsolete pairs.
        h=polys[i].LM
        C=set(basis)
        D=[]
        while C:
            j=C.pop()
            g=polys[j].LM
            L=lcm(h,g)
            masks=lambda k:quotient(L,lcm(h,polys[k].LM)) is not None
            if R.monomial_mul(h,g)==L or (not any(masks(k) for k in C) and not any(masks(k) for k in D)):
                D.append(j)
        replace=set()
        for j,k in pairs:
            u,v=polys[j].LM,polys[k].LM
            L=lcm(u,v)
            if quotient(L,h) is None or lcm(u,h)==L or lcm(v,h)==L:
                replace.add((j,k))
        for j in D:
            g=polys[j].LM
            if R.monomial_mul(h,g)!=lcm(h,g):
                replace.add((i,j))
        keep={j for j in basis if quotient(polys[j].LM,h) is None}
        keep.add(i)
        return keep,replace

    for j,v in enumerate(initial):
        i=normalize(polys[j],[(j,R.one)])
        if i is None:continue
        if polys[i].is_ground:return lift(i)
        basis,pairs=install(i)

    while pairs:
        i,j=min(pairs,key=lambda ij:order(lcm(polys[ij[0]].LM,polys[ij[1]].LM)))
        pairs.remove((i,j))
        p,q=polys[i],polys[j]
        L=lcm(p.LM,q.LM)
        u,v=quotient(L,p.LM),quotient(L,q.LM)
        left=R.from_dict({u:R.domain.one})
        right=R.from_dict({v:-R.domain.one})
        a=normalize(p.mul_monom(u)-q.mul_monom(v),[(i,left),(j,right)])
        if a is None:continue
        if polys[a].is_ground:return lift(a)
        basis,pairs=install(a)
    return None
