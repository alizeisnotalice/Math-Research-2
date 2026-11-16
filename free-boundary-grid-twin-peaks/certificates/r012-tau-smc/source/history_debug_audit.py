#!/usr/bin/env python3
"""Save full small-particle histories and independently check the two peaks.

The large GPU State keeps only a visited mask and roots.  This CPU-only audit
stores actual descending vertex orders, follows the same integer target and
fresh child-key choices, and checks the final order without consulting masks.
It is a transition audit, not an n=8 estimator or an independence claim.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

M = (1 << 64) - 1
G = 0x9E3779B97F4A7C15

def mix(x):
    z = (x + G) & M; z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M
    return z ^ (z >> 31)

def key(seed, rank, child, domain):
    return (seed ^ (rank << 32) ^ ((child * G) & M) ^ ((domain * 0xD1B54A32D192ED03) & M)) & M

def below(seed, bound):
    threshold = ((-bound) & M) % bound; state = seed
    while True:
        state = (state + G) & M; x = mix((state - G) & M)
        if x >= threshold: return x % bound

def nbrs(v, n):
    x, y = v % n, v // n
    if x: yield v - 1
    if x + 1 < n: yield v + 1
    if y: yield v - n
    if y + 1 < n: yield v + n

def eligible(path, n, rank, tau):
    seen = set(path); ans = []
    for v in range(n * n):
        if v in seen: continue
        birth = not any(w in seen for w in nbrs(v, n))
        if rank == 1 or (rank == tau and birth) or (rank != 1 and rank != tau and not birth): ans.append(v)
    return ans

def make_path(n, tau, P, seed):
    states = [()] * P; N = n * n
    for rank in range(1, N + 1):
        a = [eligible(s, n, rank, tau) for s in states]; weights = list(map(len, a)); total = sum(weights)
        if total == 0: return None
        if rank == N:
            assert all(w == 1 for w in weights)
            return list(states[0] + (a[0][0],))
        ends = []; s = 0
        for w in weights: s += w; ends.append(s)
        nxt = []
        for child in range(P):
            u = below(key(seed, rank, child, 1), total)
            parent = next(i for i, end in enumerate(ends) if u < end)
            pick = below(key(seed, rank, child, 2), weights[parent])
            nxt.append(states[parent] + (a[parent][pick],))
        states = nxt
    raise AssertionError

def independently_check(path, n, tau):
    assert len(path) == n*n == len(set(path))
    when = {v:r for r,v in enumerate(path)}
    peaks = [v for v in path if all(when[v] < when[w] for w in nbrs(v,n))]
    seen = set(); birth_ranks = []
    for r,v in enumerate(path,1):
        if not any(w in seen for w in nbrs(v,n)): birth_ranks.append(r)
        seen.add(v)
    assert birth_ranks == [1,tau] and peaks == [path[0],path[tau-1]]
    dx,dy = path[0]%n-path[tau-1]%n,path[0]//n-path[tau-1]//n
    return {"K":len(peaks),"peak_vertices":peaks,"birth_ranks":birth_ranks,"root_d2":dx*dx+dy*dy}

def main():
    n,P,R = 4,64,2; samples=[]; no_complete_history=[]
    for tau in range(2,n*n+1):
        successes=0
        for attempt in range(128):
            seed=0xA022000000000000+n*1_000_000+tau*10_000+attempt
            path=make_path(n,tau,P,seed)
            if path is None: continue
            samples.append({"tau":tau,"replicate":successes,"attempt":attempt,"seed_hex":f"{seed:016x}","descending_vertices":path,**independently_check(path,n,tau)})
            successes += 1
            if successes == R: break
        if successes < R: no_complete_history.append({"tau":tau,"successes":successes,"attempts":128})
    Path(sys.argv[1]).write_text(json.dumps({"kind":"cpu_history_transition_audit","n":n,"particles":P,"requested_replicates_per_tau":R,"samples":samples,"no_complete_history_within_attempt_budget":no_complete_history},indent=2)+"\n")
    print(f"wrote {len(samples)} complete labelled paths; every saved path independently has K=2; incomplete strata={no_complete_history}")

if __name__ == "__main__": main()
