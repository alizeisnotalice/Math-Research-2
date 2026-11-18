"""Construct both physical prefix lengths, then verify the exact tail ledger."""
from collections import Counter
import json

def check(t):
    k=t+2*(t.bit_length()-1)
    # Encode addresses by the first t bits. Each large cell uses two prefixes.
    cells=[(t-1,i) for i in range(t-1)]
    cells += [(t,i+2*(t-1)) for i in range(t)]
    assert cells[-1][1]+1<1<<t
    counts=Counter()
    ends=[]
    for depth,prefix in cells:
        lo=prefix<<(t-depth)
        hi=(prefix+1)<<(t-depth)
        ends.append((lo,hi))
        assert (1<<(k+1-depth))%4==0
        counts[depth]+=1
    assert all(end<=start for (_,end),(start,_) in zip(sorted(ends),sorted(ends)[1:]))
    q=sum(count*(1<<depth) for depth,count in counts.items())
    assert q==(3*t-1)<<(t-1) and q>3*(1<<t)-4
    return dict(t=t,k=k,r=sum(counts.values()),q=q,s=q+(1<<(k+1)),
                layer_bound=3*(1<<t)-4)
print(json.dumps([check(t) for t in range(4,65)],indent=2))
