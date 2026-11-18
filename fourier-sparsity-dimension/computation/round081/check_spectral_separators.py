"""Run python3 -B check_spectral_separators.py. Integer Walsh and GF(2) operations."""
from collections import Counter
from functools import lru_cache
import json
import platform


def walsh(values):
    out = list(values)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2*width):
            for j in range(start, start+width):
                a, b = out[j], out[j+width]
                out[j], out[j+width] = a+b, a-b
        width *= 2
    return out


def extend(H, v):
    return H | frozenset(w ^ v for w in H)


@lru_cache(None)
def span(labels):
    H = frozenset({0})
    for v in labels:
        if v not in H:
            H = extend(H, v)
    return H


def dim(H):
    return len(H).bit_length()-1


@lru_cache(None)
def all_splits(E):
    result = []
    for bits in range(1, 2**(len(E)-1)):
        B = tuple(x for i, x in enumerate(E) if (bits >> i) & 1)
        A = tuple(x for x in E if x not in B)
        if min(len(A), len(B)) < 2:
            continue
        J = span(A) & span(B)
        assert len(J) >= 2
        if len(J) == 2:
            q = next(iter(J-{0}))
            result.append((q, {min(v,v^q) for v in A}-{0}, {min(v,v^q) for v in B}-{0}))
    return result


def check_ordinary(n):
    N = 2**n
    counts = Counter()
    supports = set()
    for table in range(2**N):
        w = walsh([1-2*((table >> x)&1) for x in range(N)])
        E = tuple(e for e in range(1,N) if w[e])
        if len(E) < 4:
            continue
        d = dim(span(E))
        splits = all_splits(E)
        if splits:
            supports.add(E)
            counts["tables_with_2_separation"] += 1
        for q, A, B in splits:
            child = [{v for v in range(1,N) if v < (v^q) and w[v]+s*w[v^q]}
                     for s in (1,-1)]
            assert child == [A,B] or child == [B,A]
            s0, s1 = map(len,child)
            d0, d1 = (dim(span(tuple(sorted(c)))) for c in child)
            delta = int(w[q] != 0)
            assert len(E) == 2*(s0+s1)+delta and d == d0+d1+1
            child_full = [set(c) for c in child]
            for i,sgn in enumerate((1,-1)):
                if w[0]+sgn*w[q]:
                    child_full[i].add(0)
            ds = [dim(span(tuple(sorted(2*v+1 for v in c))))+1 for c in child_full]
            full = [v for v,a in enumerate(w) if a]
            Dp = dim(span(tuple(sorted(2*v+1 for v in full))))+1
            assert any(Dp*Dp*len(c) <= z*z*len(full) for c,z in zip(child_full,ds))
            counts["affine_potential_inequalities"] += 1
            for C in (1, 2, 7):
                D = lambda s,r: C*C*s-(r+1)**2
                assert D(len(E),d) == 2*(D(s0,d0)+D(s1,d1))+(d0-d1)**2+C*C*delta
            counts["table_partition_pairs"] += 1
    return {"variables":n,"counts":dict(counts),"different_supports_with_2_separation":len(supports)}


def interface_fixture(simultaneous):
    # q,t,guard_x,a,guard_y,b,c,d -- t is main selector, q the common interface.
    n = 8
    N = 2**n
    qbit,tbit,xbit,abit,ybit,bbit,cbit,dbit = (1<<j for j in range(n))
    def sign(x):
        q,t,xg,a,yg,b,c,d = (int(bool(x & (1 << j))) for j in range(n))
        if simultaneous:
            exponent = (xg ^ (q & a)) if t == 0 else (yg ^ (q & b) ^ (c & d))
        else:
            exponent = (xg ^ (a & t)) if q == 0 else (yg ^ (c & d))
        return 1-2*exponent
    w = walsh(map(sign,range(N)))
    E = {N+v for v,a in enumerate(w) if a}
    if simultaneous:
        H = span(tuple(sorted(v for v in E if v & xbit and not(v & tbit))))
    else:
        H = frozenset({0})
    A = {v for v in E-H if v & xbit}
    B = E-H-A
    J = span(tuple(sorted(A))) & span(tuple(sorted(B)))
    assert J == {0,qbit}
    outer_core = extend(H,qbit)
    stats = []
    for u in (0,1):
        def alive(C):
            surviving = set()
            for v in C-outer_core:
                tag = min(v,v^qbit)
                if w[tag-N]+(-1)**u*w[(tag^qbit)-N]:
                    surviving.add(tag)
            return surviving
        branches = [alive(A),alive(B)]
        both = all(branches)
        assert both == simultaneous
        projected = [{min(v^h for h in outer_core) for v in side} for side in branches]
        union = span(tuple(sorted(set.union(*projected))))
        assert dim(union) >= 2
        if both:
            # Enumerate the actual homogeneous truth function, group H' answers,
            # and test the b=t cut at every address.
            groups = {}
            for raw in range(2*N):
                if int(bool(raw & qbit)) != u:
                    continue
                key = tuple((h & raw).bit_count()%2 for h in sorted(H))
                qanswer = int(bool(raw & tbit))
                pair = groups.setdefault(key,[set(),set()])
                pair[qanswer].add(sign(raw%N)*(-1)**int(raw>=N))
            assert all(min(map(len, pair)) == 1 for pair in groups.values())
        stats.append({"interface_answer":u, "joint":bool(both), "quotient_rank":dim(union)})
    return {"mode":"joint" if simultaneous else "exclusive", "intersection_dimension":dim(J),
            "core_dimension":dim(H),"sparsity":len(E),"answers":stats}


if __name__ == "__main__":
    print(json.dumps({"python":platform.python_version(),"seed":None,"tolerance":0,
                      "exhaustion":[check_ordinary(n) for n in (3,4)],
                      "interfaces":[interface_fixture(flag) for flag in (False,True)]},indent=2))
