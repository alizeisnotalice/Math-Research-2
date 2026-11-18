"""Audit the repair by independent edge and list enumeration."""
from collections import Counter, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from itertools import permutations
from math import factorial
import json
from pathlib import Path
import platform
import subprocess
import sys


def graph(n):
    def neighbors(v):
        x, y = v % n, v // n
        return tuple(a + n*b for a, b in ((x-1,y), (x+1,y), (x,y-1), (x,y+1))
                     if 0 <= a < n and 0 <= b < n)
    return [neighbors(v) for v in range(n*n)]


def distances(n, u, v):
    return abs(u % n - v % n) + abs(u // n - v // n)


def peaks(adj, position):
    hidden = set()
    for v, neighbors in enumerate(adj):
        for u in neighbors:
            if v < u:
                hidden.add(u if position[v] < position[u] else v)
    return set(range(len(adj))) - hidden


def component(adj, seed, allowed):
    seen = {seed}
    front = [seed]
    while front:
        for v in adj[front.pop()]:
            if v in allowed and v not in seen:
                seen.add(v)
                front.append(v)
    return seen


def merge(adj, word, q):
    seen = set()
    p = word[0]
    for v in word:
        neighbors = set(adj[v]) & seen
        if q in seen and len(neighbors) >= 2:
            low = component(adj, q, seen)
            if p not in low and neighbors & low and neighbors - low:
                return v, low
        seen.add(v)
    raise AssertionError('no merge')


def shell(adj, B):
    return {w for v in B for w in adj[v] if w not in B}


def possible_components(adj, output, w, size):
    W = output[w]
    candidates = set()
    below = set()
    pending = {u for u in adj[w] if output[u] > W}
    for v in sorted((v for v in range(len(adj)) if output[v] > W),
                    key=lambda v: output[v]):
        below.add(v)
        for u in pending & below:
            c = component(adj, u, below)
            if len(c) >= size:
                pending.remove(u)
                if len(c) == size:
                    candidates.add(frozenset(c))
        if not pending:
            break
    return candidates


def branch_order(adj, B, w):
    seen = set()
    front = deque([w])
    order = []
    while front:
        for v in adj[front.popleft()]:
            if v in B and v not in seen:
                seen.add(v)
                order.append(v)
                front.append(v)
    assert seen == B
    return order


def repair_discrete(adj, word, w, B):
    bare = [v for v in word if v not in B]
    place = bare.index(w) + 1
    order = branch_order(adj, B, w)
    return bare[:place] + order + bare[place:]


def n3_exhaustive():
    n = 3
    adj = graph(n)
    bins = Counter()
    geometry = set()
    for word in permutations(range(n*n)):
        pos = {v:i+1 for i,v in enumerate(word)}
        original = peaks(adj, pos)
        if len(original) != 2:
            continue
        p = word[0]
        q, = original - {p}
        w, B = merge(adj, word, q)
        b = len(B)
        i = pos[w]
        h = min((pos[v] for v in shell(adj, B) - {w}), default=n*n+1) - i
        assert h >= 1
        bins[(i, b, h, distances(n, p, q))] += 1
        repaired = repair_discrete(adj, word, w, B)
        out = {v:j+1 for j,v in enumerate(repaired)}
        assert peaks(adj, out) == {p} and out[w] == i-b
        key = (tuple(repaired), w, b, frozenset(B))
        if key not in geometry:
            candidates = possible_components(adj, out, w, b)
            assert frozenset(B) in candidates and len(candidates) <= len(adj[w])
            geometry.add(key)
    assert sum(bins.values()) == 143808
    return bins, len(geometry)


def synthetic(n, b, fraction, move_side):
    adj = graph(n)
    N = n*n
    p = n*(n//2)+1
    w = p+1
    B = set(range(p+2, p+2+b))
    q = max(B)
    banned = B | shell(adj, B)
    chosen = [p, *sorted(B, reverse=True)]
    seen = set(chosen)
    limit = max(len(chosen), int(fraction*N))
    front = deque([p])
    while len(chosen) < limit:
        for a in adj[front.popleft()]:
            if a not in seen and a not in banned:
                seen.add(a)
                chosen.append(a)
                front.append(a)
                if len(chosen) == limit:
                    break
    chosen.append(w)
    seen.add(w)
    i = len(chosen)
    front = deque(chosen)
    while front:
        for a in adj[front.popleft()]:
            if a not in seen:
                seen.add(a)
                chosen.append(a)
                front.append(a)
    assert len(chosen) == N
    pos = {v:j*(b+2) for j,v in enumerate(chosen, 1)}
    assert peaks(adj, pos) == {p, q}
    assert merge(adj, chosen, q) == (w, B)
    M = min((pos[v] for v in shell(adj, B) - {w}), default=(N+1)*(b+2))
    W = pos[w]
    output = pos.copy()
    for j, v in enumerate(branch_order(adj, B, w), 1):
        output[v] = W+j
    assert max(output[v] for v in B) < M
    far = [v for v in range(N) if v % n >= 3*n//4 and (v%n+2*(v//n))%5 == 0]
    viable = []
    for v in far:
        largest = min(pos[u] for u in adj[v])
        if move_side == 1 and pos[v] > W and W > pos[p]+2:
            # Promote a previously below-W center to above W.
            value = min(largest-1, W-1)
        elif move_side == 0 and pos[v] > W and W+1 < largest:
            value = W+1
        else:
            continue
        # Use odd/even doubled units so the moved rank does not tie old labels.
        if value > pos[p] and distances(n, p, v) > n/8:
            viable.append((v, value))
    v, value = viable[0]
    doubled = {u:2*x for u,x in output.items()}
    doubled[v] = 2*value+1
    assert peaks(adj, doubled) == {p, v}
    assert min(doubled, key=doubled.get) == p
    rank_w = 1+sum(x < doubled[w] for x in doubled.values())
    assert rank_w == i - b + move_side
    candidates = possible_components(adj, doubled, w, b)
    assert frozenset(B) in candidates and len(candidates) <= 4
    return dict(n=n, b=b, rank_w_in=i, rank_w_out=rank_w, move_side=move_side,
                recovered_candidates=len(candidates))


def execute(binary, n, workers):
    def task(k):
        with subprocess.Popen([str(binary), str(n), str(k), str(workers)],
                              text=True, stdout=subprocess.PIPE) as process:
            out, _ = process.communicate()
            assert process.returncode == 0
        counts = Counter()
        for line in out.splitlines():
            i, b, h, d, count = map(int, line.split())
            counts[i, b, h, d] += count
        return counts

    with ThreadPoolExecutor(max_workers=workers) as pool:
        counts = sum(pool.map(task, range(workers)), Counter())
    return counts


def present(counts, n):
    by_b = defaultdict(list)
    for (i,b,h,d),v in sorted(counts.items()):
        by_b[b].append(((i,h,d), v))
    out = []
    for b, cells in sorted(by_b.items()):
        total = sum(v for _,v in cells)
        moment = sum(
            v * Fraction(factorial(h+b-1)//factorial(h-1),
                         factorial(i-1)//factorial(i-b-1))
            for (i,h,d),v in cells
        ) / total
        out.append(dict(b=b, words=total, mean_rank_ratio=str(moment)))
    return dict(n=n, A2=sum(counts.values()), by_b=out,
                histogram=[[*key, v] for key,v in sorted(counts.items())])


def main():
    b = Path(sys.argv[1]).resolve()
    data3, recovered = n3_exhaustive()
    assert execute(b, 3, 4) == data3
    data4 = execute(b, 4, 12)
    assert sum(data4.values()) == 779299174080
    synthetic_notes = [synthetic(n, k, fraction, side)
                       for n, k, fraction, side in ((64,2,0,1), (80,3,0.2,0),
                                                     (96,5,0.45,1), (128,5,0.75,1))]
    assert {1,0} <= {x['move_side'] for x in synthetic_notes}
    print(json.dumps(dict(python=platform.python_version(), n3_permutations=factorial(9),
                          checked_outputs=recovered, synthetic=synthetic_notes,
                          exact=[present(data3,3),present(data4,4)]), indent=2))


if __name__ == '__main__':
    main()
