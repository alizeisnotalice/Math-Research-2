#!/usr/bin/env python3
"""Support propagation for any degree pair of a nonnegative binary product."""

from itertools import product


def closure(a,b,positive_a,positive_b,allowed_a,allowed_b,ones,zeros,forbidden=(),required=()):
    """Return necessary support envelopes or certify the branch impossible.

    Each pair in forbidden must multiply to zero. Each family in required
    must include at least one positive product. Indices start at degree zero.
    """
    whole=set(range(a+b+1))
    pa,pb,sa,sb,one,zero=map(set,(positive_a,positive_b,allowed_a,allowed_b,ones,zeros))
    if not sa<=set(range(a+1)) or not sb<=set(range(b+1)) or not (one|zero)<=whole:
        raise ValueError("indices outside their declared degrees")
    ban=set(forbidden)
    requirements=tuple(map(set,required))
    while True:
        sizes=tuple(map(len,(pa,pb,sa,sb,one,zero)))
        active_ban=ban|{(i,j) for i,j in product(sa,sb) if i+j in zero}
        sa.difference_update(i for i,j in active_ban if j in pb)
        sb.difference_update(j for i,j in active_ban if i in pa)
        one.update(i+j for i,j in product(pa,pb))
        possible=set(product(sa,sb))-active_ban
        zero.update(whole-{i+j for i,j in possible})
        if pa-sa or pb-sb or one&zero:
            return None
        # A required positive product whose row is known zero is excluded here.
        constraints=(*requirements,*({(i,k-i) for i in sa if k-i in sb} for k in one))
        for allowed in constraints:
            known=possible&allowed
            if not known:
                return None
            if len(known)==1:
                i,j=next(iter(known))
                pa.add(i)
                pb.add(j)
        if sizes==tuple(map(len,(pa,pb,sa,sb,one,zero))):
            return (pa,pb,sa,sb,one,zero)


def saturated_pairs(a,b,saturated):
    """Translate products already equal to one into explicitly forbidden pairs."""
    chosen={}
    for i,j in saturated:
        if not 0<=i<=a or not 0<=j<=b:
            raise ValueError("saturated pair outside degree bounds")
        if i+j in chosen and chosen[i+j]!=(i,j):
            raise ValueError("two unit products in one row")
        chosen[i+j]=(i,j)
    return {(i,j) for i,j in product(range(a+1),range(b+1))
            if i+j in chosen and chosen[i+j]!=(i,j)}
