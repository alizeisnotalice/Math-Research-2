"""Check the common-basis algorithm against exhaustive small bases and actual Walsh bands."""
import itertools
import json
import platform
import random
from circuit_fields import basis


def image(v, rows):
    result = 0
    for i, r in enumerate(rows):
        if (v >> i) & 1:
            result ^= r
    return result


def all_matchings(points):
    if not points:
        yield ()
        return
    a, *rest = points
    yield from all_matchings(rest)
    for j, b in enumerate(rest):
        for tail in all_matchings(rest[:j]+rest[j+1:]):
            yield ((a,b), *tail)


def check_rows():
    rng = random.Random(8501)
    all_bases = [t for t in itertools.combinations(range(1,8),3) if len(basis(t))==3]
    count = 0
    for _ in range(512):
        maps = [[rng.randrange(8) for _ in range(3)] for _ in range(3)]
        weights = [rng.randrange(1,6) for _ in maps]
        cost = {v: sum(w for m,w in zip(maps,weights) if image(v,m)) for v in range(1,8)}
        lowest = sum(cost[v] for v in basis(sorted(cost,key=cost.get)))
        assert lowest == min(sum(cost[v] for v in b) for b in all_bases)
        compatible = lowest == sum(w*len(basis(m)) for w,m in zip(weights,maps))
        count += compatible
    return count


def matching_bands(q):
    masks = list(all_matchings(list(range(q))))
    assert len(set(masks))==len(masks)
    hist = {}
    for edges in masks:
        h = len(edges)
        hist[h] = hist.get(h,0)+1
        if q <= 8:
            vals = [(-1)**(sum(((u>>i)&1)*((u>>j)&1) for i,j in edges)%2)
                    for u in range(1<<q)]
            for b in range(q):
                step = 1<<b
                for j in range(0,len(vals),2*step):
                    for i in range(j,j+step):
                        x,y = vals[i],vals[i+step]
                        vals[i],vals[i+step] = x+y,x-y
            assert sum(w!=0 for w in vals)==4**h
            assert {abs(w) for w in vals if w}=={1<<(q-h)}
    slots = 1 << (len(masks)-1).bit_length()
    return dict(q=q, matching_count=len(masks), rank_pairs=hist, slots=slots,
                sparsity=slots*(sum(4**h*n for h,n in hist.items())+slots-len(masks)))


if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(), seed=8501, row_trials=512,
                          exact_common_basis=check_rows(),
                          matchings=[matching_bands(q) for q in (4,6,8,10)]),indent=2))
