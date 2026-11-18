"""Exact checks for the R72 simultaneous integral basis and conditional thinning."""
import collections
import itertools
import json
import math
import platform
import random
from fractions import Fraction as F
from pathlib import Path


def span(gens):
    v = {0}
    for u in gens:
        v |= {x ^ u for x in v}
    return frozenset(v)


def dot(x, v):
    return (x & v).bit_count() % 2


def inside(x, anchor, W):
    return all(dot(x ^ anchor, v) == 0 for v in W)


def main():
    X = list(range(8))
    spaces = sorted({span(J) for size in range(4)
                     for J in itertools.combinations(X[1:], size)}, key=sorted)
    flats = []
    for W in spaces:
        flats.extend((min(cs), W) for cs in sorted(
            {tuple(x for x in X if inside(x, y, W)) for y in X}))
    rng = random.Random(72026)
    M = 1000
    steps = 0
    max_cost_ratio = F(0)
    for _ in range(M):
        original = rng.sample(flats, rng.randint(1, 15))
        n = len(original)
        cur = dict(enumerate(W for _, W in original))
        weights = [rng.randint(1, 4) for _ in original]
        signed = [rng.choice((-2,-1,1,2)) for _ in original]
        mass = [F(1,len(W)) for _, W in original]
        d = [0] * n
        l = 0
        whole = span(v for W in cur.values() for v in W)
        while len(span(v for W in cur.values() for v in W)) > 1:
            total = span(v for W in cur.values() for v in W)
            order = sorted(cur)
            rng.shuffle(order)
            for i in order:
                other = span(v for j,W in cur.items() if i != j for v in W)
                if total == other:
                    del cur[i]
            new = {}
            for i, W in cur.items():
                other = span(v for j,V in cur.items() if i != j for v in V)
                U = W & other
                for v in sorted(W):
                    candidate = span([*U, v])
                    if len(candidate) < len(W):
                        U = candidate
                assert 2*len(U) == len(W)
                new[i] = U
                d[i] += 1
            cur = new
            l += 1
            steps += 1
            a = F(1,1<<l)
            residual = span(v for W in cur.values() for v in W)
            by_hash = collections.defaultdict(list)
            for x in X:
                by_hash[tuple(dot(x,v) for v in sorted(residual))].append(x)
            for cell in by_hash.values():
                count = collections.Counter(tuple(
                    i for i in cur if inside(x, *original[i])) for x in cell)
                active = [i for i,W in cur.items() if inside(cell[0], original[i][0], W)]
                for size in range(len(active)+1):
                    for S in itertools.combinations(active,size):
                        assert F(count[S],len(cell)) == a**size*(1-a)**(len(active)-size)
            for sign in (1,-1):
                num, predicted = F(0), F(0)
                for x in X:
                    num += F(2)**sum(sign*signed[i] for i in cur if inside(x,*original[i]))
                    predicted += math.prod(1+a*(F(2)**(sign*signed[i])-1)
                        for i,W in cur.items() if inside(x,original[i][0],W))
                assert num == predicted
        assert 1 << sum(d) == len(whole)
        depth = max(sum(weights[i] for i in range(n) if inside(x,*original[i])) for x in X)
        mu = sum(w*p for w,p in zip(weights,mass))
        P = sum(w*p*(1+di)**2 for w,p,di in zip(weights,mass,d))
        assert P <= mu+5*depth
        max_cost_ratio = max(max_cost_ratio,P/depth)
    print(json.dumps(dict(python=platform.python_version(), seed=72026,
          dimensions=3, flats=len(flats), families=M, verified_layers=steps,
          maximum_tested_weighted_ratio=str(max_cost_ratio),
          lower_tree_b10_k1=str(F(1,2)*(3-F(1,1024))**2),
          poisson_boundary_costs={c:(c+1)**2 for c in (2,4,8,16,32)},
          status="exact conditional cell laws, signed MGFs, rank and weighted cost verified"),
          indent=2))


if __name__ == "__main__":
    main()
