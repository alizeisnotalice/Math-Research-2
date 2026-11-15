"""Exact labelled-grid counts; Python standard library, no floating arithmetic in counts."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time


def vertices(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask ^= b


def grid(rows, cols):
    xy = [(x, y) for y in range(rows) for x in range(cols)]
    adj = [sum(1 << j for j, q in enumerate(xy)
               if abs(p[0]-q[0])+abs(p[1]-q[1]) == 1) for p in xy]
    return xy, adj


def peaks(labels, adj, region=None):
    if region is None:
        region = (1 << len(labels))-1
    return [v for v in vertices(region)
            if all(labels[v] > labels[u] for u in vertices(adj[v] & region))]


def component(mask, adj):
    seen = 0
    todo = mask & -mask
    while todo:
        b = todo & -todo
        todo ^= b
        v = b.bit_length()-1
        seen |= b
        todo |= adj[v] & mask & ~seen
    return seen


def tables(rows, cols, root_vectors=True):
    xy, adj = grid(rows, cols)
    N = len(adj)
    full = (1 << N)-1
    one, two, h = ([0]*(full+1) for _ in range(3))
    roots = [None]*(full+1) if root_vectors else None
    for S in range(1, full+1):
        if S & (S-1) == 0:
            one[S] = 1
            if roots is not None:
                r = [0]*N
                r[S.bit_length()-1] = 1
                roots[S] = r
            continue
        r = [0]*N if roots is not None else None
        for v in vertices(S):
            T = S ^ (1 << v)
            if adj[v] & T:
                one[S] += one[T]
                two[S] += two[T]
                if r is not None and one[T]:
                    for u in vertices(T):
                        r[u] += roots[T][u]
            else:
                two[S] += one[T]
        if roots is not None and one[S]:
            assert sum(r) == one[S]
            roots[S] = r
    h[full] = 1
    for S in range(full-1, 0, -1):
        h[S] = sum(h[S | (1 << v)] for v in vertices(full ^ S) if adj[v] & S)
    assert sum(h[1 << v] for v in range(N)) == one[full]
    return xy, adj, one, two, h, roots


def cut_regions(n, axis, j):
    xy, _ = grid(n, n)
    regions = [sum(1 << v for v, p in enumerate(xy) if p[axis] < j),
               sum(1 << v for v, p in enumerate(xy) if p[axis] >= j)]
    allowed = [sum(1 << v for v, p in enumerate(xy) if p[axis] <= (j-1)//2),
               sum(1 << v for v, p in enumerate(xy) if p[axis] >= j+(n-j)//2)]
    return list(zip(regions, allowed))


def event_count(n, cuts):
    """All regional restrictions are checked during descending revelation."""
    xy, adj = grid(n, n)
    full = (1 << (n*n))-1
    regions = [pair for axis, j in cuts for pair in cut_regions(n, axis, j)]
    constraints = [[(reg, allow, adj[v] & reg) for reg, allow in regions if reg >> v & 1]
                   for v in range(n*n)]
    counts = [0]*(full+1)
    counts[0] = 1
    for S in range(full):
        if not counts[S]:
            continue
        for v in vertices(full ^ S):
            b = 1 << v
            if all((bool(nb & S) if reg & S else bool(allow & b))
                   for reg, allow, nb in constraints[v]):
                counts[S | b] += counts[S]
    return counts[full]


def first_merge(order, adj):
    S, births = 0, 0
    for rank, v in enumerate(order, 1):
        if not adj[v] & S:
            births += 1
        S |= 1 << v
        if births == 2 and component(S, adj) == S:
            return rank
    return None


def brute(n):
    xy, adj = grid(n, n)
    full = (1 << (n*n))-1
    counts, dist, ranks = Counter(), Counter(), Counter()
    sq = 0
    for labels in itertools.permutations(range(1, n*n+1)):
        pp = peaks(labels, adj)
        counts[len(pp)] += 1
        if len(pp) != 2:
            continue
        p, q = (xy[v] for v in pp)
        dist[abs(p[0]-q[0])+abs(p[1]-q[1])] += 1
        sq += (p[0]-q[0])**2+(p[1]-q[1])**2
        order = sorted(range(n*n), key=lambda v: -labels[v])
        ranks[first_merge(order, adj)] += 1
    assert sum(counts.values()) == math.factorial(n*n)
    return dict(counts), dict(dist), dict(ranks), sq


def solve(n, verify_brute):
    start = time.perf_counter()
    xy, adj, one, two, h, roots = tables(n, n)
    N = n*n
    full = (1 << N)-1
    hist, rank_hist, branch_hist = Counter(), Counter(), Counter()
    sq = total = meanpart = varpart = 0
    # meanpart and varpart are exact rational weighted sums.
    meanpart = Fraction(0)
    varpart = Fraction(0)
    for S in range(1, full+1):
        if not one[S] or S.bit_count() < 3:
            continue
        for w in vertices(S):
            T = S ^ (1 << w)
            A = component(T, adj)
            B = T ^ A
            if not B or component(B, adj) != B:
                continue
            factor = h[S]*math.comb(T.bit_count(), A.bit_count())
            weight = factor*one[A]*one[B]
            total += weight
            rank_hist[S.bit_count()] += weight
            branch_hist[min(A.bit_count(), B.bit_count())] += weight
            sxA = sum(xy[u][0]*roots[A][u] for u in vertices(A))
            syA = sum(xy[u][1]*roots[A][u] for u in vertices(A))
            sxB = sum(xy[v][0]*roots[B][v] for v in vertices(B))
            syB = sum(xy[v][1]*roots[B][v] for v in vertices(B))
            md = Fraction((sxA*one[B]-sxB*one[A])**2 +
                          (syA*one[B]-syB*one[A])**2, (one[A]*one[B])**2)
            meanpart += weight*md
            local_sq = 0
            for u in vertices(A):
                for v in vertices(B):
                    c = factor*roots[A][u]*roots[B][v]
                    dx, dy = xy[u][0]-xy[v][0], xy[u][1]-xy[v][1]
                    hist[abs(dx)+abs(dy)] += c
                    local_sq += c*(dx*dx+dy*dy)
            sq += local_sq
            varpart += local_sq-weight*md
    assert total == two[full] == sum(hist.values())
    assert meanpart+varpart == sq
    if verify_brute:
        bc, bh, br, bs = brute(n)
        assert bc.get(1, 0) == one[full] and bc.get(2, 0) == two[full]
        assert bh == dict(hist) and br == dict(rank_hist) and bs == sq
    cuts = [(axis, j) for axis in (0, 1) for j in range(2, n-1)]
    singles, intersections = {}, {}
    for axis, j in cuts:
        lhs, rhs = cut_regions(n, axis, j)
        L = sum(roots[lhs[0]][u] for u in vertices(lhs[1]))
        R = sum(roots[rhs[0]][u] for u in vertices(rhs[1]))
        exact = math.comb(N, n*j)*L*R
        check = event_count(n, [(axis, j)])
        assert exact == check
        singles[f'{axis}:{j}'] = {'count':exact, 'ratio':str(Fraction(exact, total)),
            'left_single_peak_count':one[lhs[0]], 'right_single_peak_count':one[rhs[0]],
            'quarter_bound':str(Fraction(math.comb(N,n*j)*one[lhs[0]]*one[rhs[0]],4*total))}
    for c1, c2 in itertools.combinations(cuts, 2):
        intersections[f'{c1}|{c2}'] = event_count(n, [c1,c2])
    S1 = sum(v['count'] for v in singles.values())
    S2 = S1+2*sum(intersections.values())
    result = {'n':n, 'A1':one[full], 'A2':two[full],
        'distance_counts':dict(sorted(hist.items())), 'merge_rank_counts':dict(sorted(rank_hist.items())),
        'smaller_branch_counts':dict(sorted(branch_hist.items())),
        'squared_euclidean_sum':sq,
        'Q':str(Fraction(sq,N*total)) if total else None,
        'Q_mean':str(meanpart/(N*total)) if total else None,
        'Q_variance':str(varpart/(N*total)) if total else None,
        'mean_merge_rank':str(Fraction(sum(k*v for k,v in rank_hist.items()),total)) if total else None,
        'cuts':singles, 'cut_pair_intersections':intersections,
        'cut_second_moment_bound':str(Fraction(S1*S1,total*S2)) if S2 else None,
        'brute_force_verified':verify_brute, 'wall_seconds':time.perf_counter()-start}
    if len(cuts)==2:
        result['cut_union_exact_ratio'] = str(Fraction(S1-sum(intersections.values()),total))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-n',type=int,default=4)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'arithmetic':'Python int and Fraction; deterministic', 'results':[]}
    for n in range(1,args.max_n+1):
        r=solve(n,n<=3)
        out['results'].append(r)
        args.output.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(r),flush=True)


if __name__=='__main__':
    main()
