"""Certified large fibres for canonical internal-geodesic promotion.

The construction is a witness against a uniform C^d d! fibre bound only; it
does not compare its family with all two-peak labelings.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def xy(v, n):
    return v % n, v // n


def vid(x, y, n):
    return y*n+x


def neighbours(v, n):
    x,y=xy(v,n)
    return [vid(X,Y,n) for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1))
            if 0<=X<n and 0<=Y<n]


def connected_order(vertices, start, n, required_prefix):
    """Grow a connected order, forcing connected required_prefix first."""
    allowed=set(vertices); order=[start]; seen={start}
    for v in required_prefix:
        assert v in allowed and v not in seen
        assert any(w in seen for w in neighbours(v,n))
        order.append(v);seen.add(v)
    while len(order)<len(allowed):
        choices=sorted(v for v in allowed-seen if any(w in seen for w in neighbours(v,n)))
        assert choices
        v=choices[0];order.append(v);seen.add(v)
    return order


def build(n):
    assert n>=4
    N=n*n;u=vid(0,0,n);v=vid(n-1,0,n);q=vid(n-1,1,n)
    internal=[vid(x,0,n) for x in range(1,n-1)]
    anchors=[vid(x,1,n) for x in range(1,n-1)]
    H=set(range(N))-set(internal)
    # The displayed initial chain connects u to every anchor without v or q.
    initial=[vid(0,1,n)]+anchors
    allowed=H-{v,q}
    assert set(initial)<={*allowed}
    pre=connected_order(allowed,u,n,initial)
    r=N//2
    # v is the r-th H vertex.  Its only H-neighbour q remains absent.
    assert len(pre)>=r-1 and r-1>=len(initial)+1
    prefix=pre[:r-1]
    rest=pre[r-1:]
    H_order=prefix+[v,q]+rest
    assert len(H_order)==len(H) and H_order[0]==u and H_order.index(v)==r-1
    # Every non-u H vertex except v has a higher H neighbour; v is a new birth.
    revealed=set()
    births=[]
    for w in H_order:
        if not any(z in revealed for z in neighbours(w,n)): births.append(w)
        revealed.add(w)
    assert births==[u,v]
    h=len(H_order)-(H_order.index(v)+1)
    return {'n':n,'N':N,'u':u,'v':v,'q':q,'internal':internal,'anchors':anchors,
            'H_order':H_order,'tail_length':h,'fiber_size':math.factorial(h+len(internal))//math.factorial(h)}


def labels_from_order(order):
    N=len(order); labels=[0]*N
    for rank,v in enumerate(order): labels[v]=N-rank
    return labels


def peaks(labels,n):
    return [v for v in range(n*n) if all(labels[v]>labels[w] for w in neighbours(v,n))]


def promoted(order, data):
    internal=set(data['internal']);u=data['u']
    remainder=[v for v in order if v not in internal]
    assert remainder[0]==u
    return [u]+data['internal']+remainder[1:]


def sample_orders(data, exhaustive_limit=500000):
    H=data['H_order'];p=data['internal'];i=H.index(data['v']); head=H[:i+1];tail=H[i+1:]
    total=math.factorial(len(tail)+len(p))//math.factorial(len(tail))
    if total>exhaustive_limit:
        # Deterministic boundary cases, followed by two mixed placements.
        yield head+p+tail
        yield head+tail+p
        for perm in (p, list(reversed(p))):
            positions={2*k for k in range(len(p))}
            merged=[];a=iter(tail);b=iter(perm)
            for k in range(len(tail)+len(p)):
                merged.append(next(b) if k in positions else next(a))
            yield head+merged
        return
    for positions in itertools.combinations(range(len(tail)+len(p)),len(p)):
        pos=set(positions)
        for perm in itertools.permutations(p):
            a=iter(tail);b=iter(perm)
            yield head+[next(b) if k in pos else next(a) for k in range(len(tail)+len(p))]


def audit(n):
    data=build(n);outputs=set();checked=0
    for order in sample_orders(data):
        labels=labels_from_order(order)
        assert peaks(labels,n)==[data['u'],data['v']]
        # Each path vertex is below its pre-revealed anchor in the source.
        for z,a in zip(data['internal'],data['anchors']): assert labels[a]>labels[z]
        out=promoted(order,data);outlabels=labels_from_order(out)
        assert peaks(outlabels,n)==[data['u']]
        assert out[:n-1]==[data['u']]+data['internal']
        outputs.add(tuple(out));checked+=1
    assert len(outputs)==1
    data['audited_orders']=checked;data['exhaustive']=data['fiber_size']<=500000
    data['common_output']=list(next(iter(outputs)))
    return data


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[4,5]);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out={'results':[audit(n) for n in a.n]};a.output.write_text(json.dumps(out,indent=2)+'\n')
    for r in out['results']:print(json.dumps({k:r[k] for k in ('n','N','tail_length','fiber_size','audited_orders','exhaustive')}))
if __name__=='__main__':main()
