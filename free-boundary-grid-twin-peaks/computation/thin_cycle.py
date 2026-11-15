"""Verify the 3-by-n induced-boundary counterexample in the n-by-n grid."""
from fractions import Fraction


def adj(n, p, q):
    return abs(p[0] - q[0]) + abs(p[1] - q[1]) == 1


def cycle_vertices(n):
    return {(x, y) for y in range(n) for x in range(3) if x in (0, 2) or y in (0, n - 1)}


def verify(n):
    A = cycle_vertices(n)
    degrees = {p: sum(adj(n, p, q) for q in A) for p in A}
    assert len(A) == 2 * n + 2
    assert set(degrees.values()) == {2}
    # A finite connected 2-regular graph is a cycle; explicitly check connectivity.
    reached = {next(iter(A))}
    while True:
        new = reached | {q for p in reached for q in A if adj(n, p, q)}
        if new == reached:
            break
        reached = new
    assert reached == A
    ybar = Fraction(n - 1, 2)
    vy = sum((Fraction(y) - ybar) ** 2 for _, y in A) / len(A)
    claimed = Fraction((n - 1) * (n * n + 4 * n - 3), 12 * (n + 1))
    assert vy == claimed
    r = n // 2
    w, B = (3, r), {(4, r)}
    assert n >= 5 and w not in A and B.isdisjoint(A)
    assert sum(adj(n, w, q) for q in A) == 1
    assert sum(adj(n, w, q) for q in B) == 1
    assert not any(adj(n, p, q) for p in A for q in B)
    return {"n": n, "A_size": len(A), "root_y_variance": str(vy),
            "root_y_variance_float": float(vy), "w": w, "B": sorted(B)}



def weight_n5():
    from functools import lru_cache
    from math import comb
    from exact import grid, vertices
    n=5
    xy,nb=grid(n,n)
    A=sum(1<<v for v,q in enumerate(xy) if q in cycle_vertices(n))
    B=1<<(2*n+4);w=2*n+3;S=A|B|(1<<w);full=(1<<25)-1
    @lru_cache(None)
    def one(T):
        if T and T&(T-1)==0:return 1
        return sum(one(T^(1<<v)) for v in vertices(T) if nb[v]&(T^(1<<v)))
    @lru_cache(None)
    def complete(T):
        if T==full:return 1
        return sum(complete(T|(1<<v)) for v in vertices(full^T) if nb[v]&T)
    a=one(A);h=complete(S);W=h*comb(S.bit_count()-1,A.bit_count())*a
    assert (a,h,W)==(12288,17740800,2833986355200)
    return {"n":n,"a_A":a,"h_S":h,"W":W,"completion_states":complete.cache_info().currsize}

if __name__=="__main__":
    import json
    print(json.dumps({"variance_checks":[verify(n) for n in (5,7,11,31)],"weight_check":weight_n5()},indent=2))
