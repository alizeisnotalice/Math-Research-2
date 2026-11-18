"""Finite-field equation construction and exact affine substitutions."""
from prefix_29_30 import slots,model
from verify_degrees_25_26 import make_ring,truncated_product

def equations(j,shape,mod):
    prime,t,weight=mod
    R,*vs=make_ring(prime,len(slots(j,shape)))
    z=vs[-1]
    E,A,B=model(R(t),j,shape,vs[:-1])
    rhs=sum(w*truncated_product(A,B,k) for k,w in enumerate(weight,j-7))
    return R,vs,[R(s) for _,s in E]+[z-rhs]


def substitution(R,f,i,L,x):
    """Return f(L) and an exact quotient of [f(x)-f(L)] by x-L."""
    total=R.zero
    quotient=R.zero
    powers=[R.one]
    for powers_xyz,c in f.items():
        e=powers_xyz[i]
        base=R.from_dict({powers_xyz[:i]+(0,)+powers_xyz[i+1:]:c})
        while len(powers)<=e:
            powers.append(powers[-1]*L)
        total+=base*powers[e]
        if e:
            quotient+=base*sum((x**(e-1-k))*powers[k] for k in range(e))
    assert f-total==quotient*(x-L)
    return total,quotient


def forward_eliminate(R,vars,E,track):
    pairs=[(g,[R.one if i==k else R.zero for k in range(len(E))] if track else None) for i,g in enumerate(E) if g]
    while True:
        unit=next((x for x in pairs if x[0].is_ground),None)
        if unit is not None:
            return [unit]
        best=None
        for k,(g,d) in enumerate(pairs):
            if any(sum(m)>1 for m in g):
                continue
            for m,v in g.items():
                if any(m[:-1]):
                    q=next(i for i,x in enumerate(m) if x)
                    rank=(len(g),k,q)
                    if best is None or rank < best[0]:
                        best=(rank,v)
        if best is None:
            return pairs
        (_,k,i),coef=best
        g,d=pairs.pop(k)
        u=R.domain.one/coef
        x=vars[i]
        L=-(g*u-x)
        factor=[f*u for f in d] if track else None
        revised=[]
        for h,provenance in pairs:
            if all(not e[i] for e in h):
                revised.append((h,provenance))
                continue
            y,q=substitution(R,h,i,L,x)
            if y:
                next_prov=[b-q*f for b,f in zip(provenance,factor,strict=True)] if track else None
                revised.append((y,next_prov))
        pairs=revised


