"""Add the disjoint single-dominant-winner quadrant family exactly."""
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import comb, factorial
from pathlib import Path

from merge_quadrants import HERE, Ranks, adjacency, decimal, sha


def read_ranks(m):
    data = json.loads((HERE/f'double-m{m}.json').read_text())
    assert data['m'] == m
    return {row['root']: [[[0,*map(int, s)] for s in r] for r in row['last_connector_rank_counts']]
            for row in data['roots']}


def brute(m, recorded):
    area = m*m
    adj = adjacency(m)
    response = {u: [[[0]*(area+1) for _ in range(m-1)] for _ in range(m-1)] for u in recorded}
    for word in permutations(range(area)):
        u = word[0]
        if u not in response:
            continue
        prefix = 1<<u
        for v in word[1:]:
            if not adj[v] & prefix:
                break
            prefix |= 1<<v
        else:
            rank = {v:j for j,v in enumerate(word,1)}
            for a in range(m-1):
                for b in range(m-1):
                    j = max(rank[(m-1)*m+a], rank[b*m+m-1])
                    response[u][a][b][j] += 1
    assert response == recorded


def triple_shuffles(area,j):
    return factorial(3*area-j)//(factorial(area-j)*factorial(area)**2)


def audit_shuffles(area):
    for j in range(1,area+1):
        # Choose the dominant block's positions and then unconditionally split
        # the remainder into the two losing blocks.
        all_count = sum(comb(2*area,area) for A in combinations(range(3*area),area) if A[j-1]==j-1)
        assert triple_shuffles(area,j) == all_count


def build_core(m, root, joint):
    area = m**2
    weights = [0] + [triple_shuffles(area,j) for j in range(1,area+1)]
    for u,bins in joint.items():
        uT = (u%m)*m+u//m
        for a in range(m-1):
            for b in range(m-1):
                assert sum(bins[a][b]) == root[u]
                assert bins[a][b] == joint[uT][b][a]
    dominant = {u:sum(root[a]*root[b]*sum(c*d for c,d in zip(hist,weights))
                      for a,row in enumerate(bins) for b,hist in enumerate(row)) for u,bins in joint.items()}
    data = []
    n = 2*m
    for u,v in combinations(range(n*n),2):
        x,y,z,w = u%n,u//n,v%n,v//n
        if x//m == z//m or y//m == w//m:
            continue
        xc,yc,zc,wc = (min(h,n-1-h) for h in (x,y,z,w))
        if max(xc,yc,zc,wc)>=m-1:
            continue
        uc,vc = yc*m+xc,wc*m+zc
        count = comb(4*area,area)*(dominant[uc]*root[vc]+root[uc]*dominant[vc])
        assert count > 0
        data.append({'peaks':[u,v], 'manhattan_distance':abs(x-z)+abs(y-w), 'new_count':str(count)})
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('axis_union', type=Path)
    parser.add_argument('output', type=Path)
    a = parser.parse_args()
    axis = json.loads(a.axis_union.read_text())
    for m in (2,3):
        brute(m,read_ranks(m))
    for area in (1,2,3,4):
        audit_shuffles(area)
    m=5
    book=Ranks(HERE/'ranks-m5.json')
    extra=build_core(m,book.root,read_ranks(m))
    bound=int(axis['denominator_A2_upper'])
    incumbent=int(axis['incumbent_count_reproduced'])
    added_by_dist=defaultdict(int)
    for row in extra:
        added_by_dist[row['manhattan_distance']]+=int(row['new_count'])
    new_thresholds=[]
    for old in axis['paired_seam_roots_also_allowed']['thresholds']:
        q=Fraction(old['strict_manhattan_threshold'])
        added=sum(c for d,c in added_by_dist.items() if d>q)
        total=int(old['union_count'])+added
        assert total<=bound
        lower=Fraction(total,bound)
        new_thresholds.append({'strict_manhattan_threshold':str(q),'new_disjoint_count':str(added),
                               'axis_union_count':old['union_count'],'combined_exact_count':str(total),
                               'conditional_lower_exact':str(lower),'conditional_lower_decimal':decimal(lower),
                               'factor_over_incumbent':decimal(Fraction(total,incumbent)),
                               'factor_over_axis_bound':decimal(Fraction(total,int(old['union_count'])))})
    output={'scope':'Exact additional class disjoint from axis union, finite bound only',
            'source_sha256':sha(HERE/'double_connectors.cpp'),'rank_ledger_sha256':sha(HERE/'double-m5.json'),
            'axis_ledger_sha256':sha(a.axis_union),'denominator_A2_upper':str(bound),
            'additional_root_pairs':len(extra),'additional_per_dist':{str(d):str(c) for d,c in sorted(added_by_dist.items())},
            'added_strata':extra,'thresholds':new_thresholds}
    with a.output.open('x') as f:
        json.dump(output,f,indent=2)
        f.write('\n')
    for row in new_thresholds:
        print(row['strict_manhattan_threshold'],row['conditional_lower_decimal'],row['factor_over_axis_bound'])

if __name__=='__main__':
    main()
