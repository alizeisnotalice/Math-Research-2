"""Audit a proposed top-two inclusion-exclusion identity on the 3x3 grid."""
from fractions import Fraction
from itertools import permutations
from math import factorial

def audit(adj):
    N = len(adj)
    full = (1 << N) - 1
    edges = [(v,u) for v in range(N) for u in range(v+1,N) if adj[v] >> u & 1]
    closed = [adj[v] | (1 << v) for v in range(N)]

    def is_independent(mask):
        return all(not (mask >> v & 1 and adj[v] & mask) for v in range(N))

    independent = sorted((s for s in range(1 << N) if is_independent(s)), key=int.bit_count)
    C = {0: factorial(N)}
    for I in independent[1:]:
        q = I
        neighbourhood = I
        numerator = 0
        while q:
            b = q & -q
            q -= b
            numerator += C[I ^ b]
            neighbourhood |= closed[b.bit_length()-1]
        C[I], rem = divmod(numerator, neighbourhood.bit_count())
        assert rem == 0

    def edge_count(mask):
        return sum((mask >> a & 1) and (mask >> b & 1) for a,b in edges)

    candidate_adj = Fraction()
    candidate_nonadj = Fraction()
    for I in independent:
        sign = -1 if I.bit_count() & 1 else 1
        Z = full
        q = I
        while q:
            b = q & -q
            q -= b
            Z &= ~closed[b.bit_length()-1]
        z = Z.bit_count()
        e = edge_count(Z)
        candidate_adj += sign * C[I] * Fraction(2*e, N*(N-1))
        candidate_nonadj += sign * C[I] * Fraction(z*(z-1)-2*e, N*(N-1))

    direct = {'K1':0, 'K2':0, 'top2_adj':0, 'top2_nonadj':0,
              'K1_top2_adj':0, 'K1_top2_nonadj':0,
              'K2_top2_adj':0, 'K2_top2_nonadj':0}
    for labels in permutations(range(1,N+1)):
        peaks = [v for v in range(N) if all(labels[v] > labels[u] for u in range(N) if adj[v] >> u & 1)]
        top = labels.index(N)
        second = labels.index(N-1)
        category = 'top2_adj' if adj[top] >> second & 1 else 'top2_nonadj'
        direct[category] += 1
        if len(peaks) == 1:
            direct['K1'] += 1; direct['K1_'+category] += 1
        if len(peaks) == 2:
            direct['K2'] += 1; direct['K2_'+category] += 1
    return {'candidate_adj': str(candidate_adj), 'candidate_nonadj': str(candidate_nonadj),
            'direct': direct, 'adj_ok':candidate_adj==direct['K1'],
            'nonadj_ok':candidate_nonadj==direct['K2_top2_nonadj']}

grid=[0]*9
for y in range(3):
    for x in range(3):
        v=y*3+x
        for X,Y in ((x-1,y),(x+1,y),(x,y-1),(x,y+1)):
            if 0<=X<3 and 0<=Y<3:grid[v]|=1<<(Y*3+X)
path=[0]*4
for v in range(3):path[v]|=1<<(v+1);path[v+1]|=1<<v
print({'grid3':audit(grid),'path4':audit(path)})
