"""python3 -B check_unrestricted.py. Exact adversarial weight checks for unrestricted.md."""
import json
import platform
import random
from collections import Counter
from itertools import combinations, product


def rank(vs):
    piv={}
    for v in vs:
        while v and v.bit_length() in piv:
            v ^= piv[v.bit_length()]
        if v:
            piv[v.bit_length()]=v
    return len(piv)


def columns(A,a,b):
    return [(A>>(a*j))&((1<<a)-1) for j in range(b)]


def rows(A,a,b):
    cs=columns(A,a,b)
    return [sum(((c>>j)&1)<<i for i,c in enumerate(cs)) for j in range(a)]


def alt(mask,n):
    out=0
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if (mask>>k)&1:
            out|=(1<<(n*i+j))|(1<<(n*j+i))
    return out


def min_cover(triple):
    # Extreme points of the pair-cover constraints w_i+w_j>=2 and w_i>=0.
    return min(sum(triple),2*(sum(triple)-max(triple)))


def examine(A,B,a,b,alternating):
    mats=(A,B,A^B)
    cs=[columns(C,a,b) for C in mats]
    ranks=list(map(rank,cs))
    t=rank([v for group in cs for v in group])
    if not t:
        return (0,0)
    if alternating:
        assert all(r%2==0 for r in ranks)
        hs=[r//2 for r in ranks]
        target=t
    else:
        hs=ranks
        target=t+rank([v for M in mats for v in rows(M,a,b)])
    f=[(1<<h)-1 for h in hs]
    den=1<<max(hs)
    g=[(1<<h)*den-(den>>h) for h in hs]
    slack_f=min_cover(f)-target
    slack_g=min_cover(g)-(target+1)*den
    assert min(slack_f,slack_g)>=0,(mats,ranks,slack_f,slack_g)
    return min(1,slack_f),min(1,slack_g)


def counts(matrices,a,b,alternating):
    records=Counter(examine(A,B,a,b,alternating) for A,B in matrices)
    return {'families':sum(records.values()),'both_rates_tight_or_zero':records[(0,0)],
            'h_tight_phi_loose':records[(0,1)],'h_loose_phi_tight':records[(1,0)],
            'both_strict':records[(1,1)]}


def main():
    checks={}
    for a,b in [(1,4),(2,2),(2,3),(3,2)]:
        checks[f'rectangular_{a}x{b}']=counts(product(range(1<<(a*b)),repeat=2),a,b,False)
    for n in (3,4):
        catalog=[alt(i,n) for i in range(1<<(n*(n-1)//2))]
        checks[f'alternating_{n}']=counts(product(catalog,repeat=2),n,n,True)
    rng=random.Random(20261118)
    for n in (5,6,8,12):
        width=n*(n-1)//2
        pairs=((alt(rng.getrandbits(width),n),alt(rng.getrandbits(width),n)) for _ in range(4096))
        checks[f'alternating_sample_{n}']=counts(pairs,n,n,True)
    print(json.dumps({'python':platform.python_version(),'tolerance':0,'seed':20261118,
                      'data':'all nonzero vectors of F2^2, all admissible nonnegative cut weights',
                      'checks':checks},indent=2))


if __name__=='__main__':
    main()
