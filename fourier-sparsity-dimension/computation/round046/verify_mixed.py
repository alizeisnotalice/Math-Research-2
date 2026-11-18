"""Deterministic exact small tests of original parity trees and their quotient rewrites."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
import platform
import random


def bit(a, b):
    return (a & b).bit_count() % 2


def spectrum(values):
    out = list(values)
    h = 1
    while h < len(out):
        for i in range(0, len(out), 2*h):
            for j in range(i, i+h):
                u, v = out[j], out[j+h]
                out[j], out[j+h] = u+v, u-v
        h *= 2
    return {i for i, a in enumerate(out) if a}


def span(vs):
    S = {0}
    for v in vs:
        S |= {a ^ v for a in S}
    return S


def make(seed):
    rng = random.Random(seed)
    def subtree(path):
        if len(path) == 4:
            return rng.getrandbits(8)
        if not path:
            v = 1
        elif len(path) == 1:
            v = (1, 2)[path[0]]
        elif len(path) == 2:
            v = (2, 2, 3, 4)[2*path[0]+path[1]]
        else:
            v = rng.choice((1, 2, 3, 4, 5, 6, 7))
        q = v | (rng.randrange(4) << 3)
        return q, subtree(path+(0,)), subtree(path+(1,))
    return subtree(())


def audit(seed):
    tree = make(seed)
    nodes = {}
    def rewrite(node, path, private, feasible):
        if isinstance(node, int):
            return ('l', node)
        q, left, right = node
        qs = list(private)
        for mask in range(1 << len(qs)):
            s, c = q, 0
            for j,(ancestor, side) in enumerate(qs):
                if (mask >> j) & 1:
                    s ^= ancestor
                    c ^= side
            if not (s & 7):
                shared = s >> 3
                nodes[path] = ('s', shared, c, qs, feasible)
                sub = [rewrite(child, path+(side,), qs,
                               {w for w in feasible if (bit(w, shared)^c)==side})
                       for side,child in enumerate((left,right))]
                return ('s', shared, c, *sub)
        nodes[path] = ('p', q & 7, q >> 3, qs, feasible)
        return ('p', q & 7, q >> 3, *(rewrite(child, path+(side,), qs+[(q,side)],feasible)
                                         for side,child in enumerate((left,right))))
    root = rewrite(tree,(),[],set(range(4)))
    def raw(node,x,w,y):
        if isinstance(node,int):
            return (-1)**(node >> (w+4*y) & 1)
        q,l,r=node
        return raw((l,r)[bit(q,x | w<<3)],x,w,y)
    def norm(node,x,w,y):
        if node[0]=='l':
            return (-1)**(node[1] >> (w+4*y) & 1)
        typ,q,c,l,r=node
        a = bit(q,x)^bit(c,w) if typ=='p' else bit(q,w)^c
        return norm((l,r)[a],x,w,y)
    def mean(node,w,y):
        if node[0]=='l':
            return (-1)**(node[1] >> (w+4*y)&1)
        typ,q,c,l,r=node
        if typ=='s':
            return mean((l,r)[bit(q,w)^c],w,y)
        return (mean(l,w,y)+mean(r,w,y))/Fraction(2)
    YX=[(x,w,y) for y in range(2) for w in range(4) for x in range(8)]
    values=[raw(tree,x,w,y) for x,w,y in YX]
    assert values==[norm(root,x,w,y) for x,w,y in YX]
    actual=spectrum(values)
    kraft=defaultdict(Fraction)
    bound=Fraction(0)
    count=Counter()
    for path,(typ,q,c,qs,ws) in nodes.items():
        if not ws: continue
        if typ=='s':
            count['converted']+=1
            continue
        ancestors=[u & 7 for u,_ in qs]
        subtree=root
        for side in path:subtree=subtree[-2+side]
        L,R=subtree[-2:]
        hv=[(mean(L,w,y)-mean(R,w,y))/Fraction(2) if w in ws else 0
            for y in range(2) for w in range(4)]
        sh=len(spectrum(hv))
        if not sh:continue
        count['marked']+=1
        E={i ^ q for i in span(ancestors)}
        choice=None
        for toggle in range(1,8):
            if not bit(toggle,q) or any(bit(toggle,o) for o in ancestors):continue
            blockers=set()
            for z,(zt,zq,_,_,zw) in nodes.items():
                if zt!='p' or not zw or z[:len(path)]==path or not bit(toggle,zq):continue
                i=next(i for i,(a,b) in enumerate(zip(z,path)) if a!=b)
                blockers.add(path[:i])
            a=sum(nodes[z][0]=='p' for z in blockers)
            shared=[nodes[z][1] for z in blockers if nodes[z][0]=='s']
            b=len(shared); r=len(span(shared)).bit_length()-1
            item=(2*a+b+r,a,b,r)
            if choice is None or item<choice:choice=item
        _,a,b,r=choice
        count['marked_with_shared_blocker']+=int(b!=0)
        assert sum((xi&7) in E for xi in actual) * (1<<(a+r)) >= (1<<len(qs))*sh
        for xi in E:kraft[xi]+=Fraction(1,1<<(a+b))
        bound+=Fraction((1<<len(qs))*sh,1<<(2*a+b+r))
    assert max(kraft.values(),default=0)<=1 and bound<=len(actual)
    return len(actual),bound,count


tot=Counter()
for seed in range(400):
    s,b,c=audit(seed)
    tot.update(c)
    tot['trees']+=1
    tot['joint_coefficients']+=64
    tot['support_total']+=s
    tot['inequality_equalities']+=int(s==b)
assert tot['marked_with_shared_blocker']>0
print(json.dumps(dict(python=platform.python_version(), results=tot),indent=2))
