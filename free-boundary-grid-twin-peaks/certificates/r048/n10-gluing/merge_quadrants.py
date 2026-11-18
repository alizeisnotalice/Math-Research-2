"""Exact one-peak quadrant gluing, with all label interleavings and axis overlap."""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal(ratio):
    with localcontext() as ctx:
        ctx.prec = 55
        return str(Decimal(ratio.numerator) / Decimal(ratio.denominator))


class Ranks:
    def __init__(self, path):
        ledger = json.loads(path.read_text())
        self.m = ledger['m']
        self.area = self.m**2
        self.root = list(map(int, ledger['rooted_counts']))
        self.rank = {row['root']: [[0, *map(int, h)] for h in row['rank_counts']]
                     for row in ledger['roots']}
        for u, ranks in self.rank.items():
            assert len(ranks) == self.m and all(len(row) == self.area+1 and sum(row) == self.root[u] for row in ranks)

    def ranks(self, x, y, connector):
        if x > (self.m-1)//2:
            x = self.m-1-x
            connector = self.m-1-connector
        return self.rank[y*self.m+x][connector]

    def pair_counts(self):
        return {(x,y):sum(self.root[c]*a*comb(2*self.area-j,self.area)
                          for c in range(self.m)
                          for j,a in enumerate(self.ranks(x,y,c)) if a)
                for x in range(self.m-1) for y in range(self.m)}


def overlap_shuffles(area):
    fact = [factorial(i) for i in range(4*area+1)]
    suffix_sum = [[0]*(area+2) for _ in range(area+2)]
    for a in range(area,-1,-1):
        for b in range(area,-1,-1):
            suffix_count = fact[4*area-a-b-1]//(fact[area-a]*fact[area-b]*fact[area-1]*fact[area])
            value = 2*comb(a+b,a)*suffix_count
            suffix_sum[a][b] = value+suffix_sum[a+1][b]+suffix_sum[a][b+1]-suffix_sum[a+1][b+1]
    return suffix_sum


def brute_overlap_shuffles(area):
    counts = [[0]*(area+2) for _ in range(area+2)]
    vertices = set(range(4*area))
    for a in combinations(range(4*area),area):
        unused = vertices.difference(a)
        for b in combinations(sorted(unused),area):
            suffix = unused.difference(b)
            first_loser = min(suffix)
            aa = sum(x < first_loser for x in a)
            bb = sum(x < first_loser for x in b)
            counts[aa][bb] += 2*comb(2*area-1,area-1)
    return [[sum(counts[a][b] for a in range(i,area+1) for b in range(j,area+1))
             for j in range(area+2)] for i in range(area+2)]


def adjacency(size):
    return [sum(1<<u for u in range(size**2) if abs(u%size-v%size)+abs(u//size-v//size)==1)
            for v in range(size**2)]


def brute_rank_audit(book):
    m = book.m
    assert m in (2,3)
    area = book.area
    adj = adjacency(m)
    root_counts = [0]*area
    rank_counts = {u:[[0]*(area+1) for _ in range(m)] for u in book.rank}
    for word in permutations(range(area)):
        revealed = 1<<word[0]
        for v in word[1:]:
            if not adj[v] & revealed:
                break
            revealed |= 1<<v
        else:
            u = word[0]
            root_counts[u] += 1
            if u in rank_counts:
                for j,v in enumerate(word,1):
                    if v//m == m-1:
                        rank_counts[u][v%m][j] += 1
    assert book.root == root_counts and book.rank == rank_counts
    return factorial(area)


def reproduce_incumbent(book, incumbent):
    area = book.area
    root_band_width = book.m//2
    aligned_pairs = sum(book.root[x]**2 for x in range(root_band_width))
    mixed = 2*comb(2*area-1,area)*aligned_pairs
    old_words = comb(4*area,2*area)*mixed**2
    assert old_words == int(incumbent['constructed_words'])
    assert len(book.root) == len(incumbent['all_root_counts']) and book.root == incumbent['all_root_counts']
    return old_words


def geometric_union(book, allow_pair_seam):
    m, area = book.m, book.area
    n = 2*m
    pair = book.pair_counts()
    shuffle_pairs = comb(n*n,2*area)
    cross_shuffle = overlap_shuffles(area)
    central_loser = book.root[0]
    records = []
    for u,v in combinations(range(n*n),2):
        x,y = u%n,u//n
        z,w = v%n,v//n
        xc,yc,zc,wc = (min(a,n-1-a) for a in (x,y,z,w))
        if not allow_pair_seam and (m-1 in (xc,yc,zc,wc)):
            continue
        vertical = x//m != z//m and xc < m-1 and zc < m-1
        horizontal = y//m != w//m and yc < m-1 and wc < m-1
        if not (vertical or horizontal):
            continue
        col = shuffle_pairs*pair[xc,yc]*pair[zc,wc] if vertical else 0
        row = shuffle_pairs*pair[yc,xc]*pair[wc,zc] if horizontal else 0
        inter = 0
        if vertical and horizontal:
            i_ranks = book.ranks(xc,yc,m-1)
            j_ranks = book.ranks(zc,wc,m-1)
            inter = central_loser**2 * sum(ai*bj*cross_shuffle[i][j] for i,ai in enumerate(i_ranks) if ai
                                           for j,bj in enumerate(j_ranks) if bj)
            assert 0 < inter <= min(row,col)
        dist = abs(x-z)+abs(y-w)
        total = col+row-inter
        assert total > 0
        records.append({'peaks':[u,v], 'manhattan_distance':dist, 'vertical_count':str(col),
                        'horizontal_count':str(row),'intersection_count':str(inter),'union_count':str(total)})
    return records


def summary(m, records, denominator, incumbent):
    n = 2*m
    by_dist = defaultdict(int)
    results = []
    for row in records:
        by_dist[row['manhattan_distance']] += int(row['union_count'])
    full_count = sum(by_dist.values())
    assert full_count <= denominator
    for threshold in (Fraction(n,3), Fraction(n,2), Fraction(6)):
        chosen = [r for r in records if r['manhattan_distance']>threshold]
        union = sum(int(r['union_count']) for r in chosen)
        v = sum(int(r['vertical_count']) for r in chosen)
        h = sum(int(r['horizontal_count']) for r in chosen)
        both = sum(int(r['intersection_count']) for r in chosen)
        assert union == v+h-both and v==h
        lower = Fraction(union,denominator)
        results.append({'strict_manhattan_threshold':str(threshold), 'root_pairs':len(chosen),
                        'vertical_count':str(v),'horizontal_count':str(h),'intersection_count':str(both),
                        'union_count':str(union),'original_conditional_lower_exact':str(lower),
                        'original_conditional_lower_decimal':decimal(lower),
                        'factor_over_incumbent':decimal(Fraction(union,incumbent))})
    return {'n':n, 'full_family_union_count':str(full_count),'per_dist_union_count':{str(d):str(c) for d,c in sorted(by_dist.items())},
            'unconditional_union_probability_decimal':decimal(Fraction(full_count,factorial(n*n))),
            'root_pairs':len(records),'thresholds':results, 'strata':records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('incumbent',type=Path)
    parser.add_argument('denominator',type=Path)
    options = parser.parse_args()
    old = json.loads(options.incumbent.read_text())
    reference = json.loads(options.denominator.read_text())
    assert sha(options.denominator) == old['from_upper_sha256']
    A2_upper = int(reference['upper_counts'][0]['A2_upper_count'])
    for i in (1,2,3):
        assert overlap_shuffles(i) == brute_overlap_shuffles(i)
    exhaustive = sum(brute_rank_audit(Ranks(HERE/f'ranks-m{m}.json')) for m in (2,3))
    book = Ranks(HERE/'ranks-m5.json')
    safe_old = Ranks(HERE/'ranks-m5-safe-only.json')
    assert all(value == book.rank[u] for u,value in safe_old.rank.items())
    incumbent = reproduce_incumbent(book,old)
    safe = geometric_union(book,allow_pair_seam=False)
    larger = geometric_union(book,allow_pair_seam=True)
    reduced = {tuple(row['peaks']):row for row in larger if all(v%10 not in (4,5) and v//10 not in (4,5) for v in row['peaks'])}
    assert {tuple(row['peaks']):row for row in safe} == reduced
    result = {'scope':'Exact finite cut-family union, conditional lower via existing upper A2 bound',
              'rank_counter_sha256':sha(HERE/'quadrant_ranks.cpp'),'m5_rank_ledger_sha256':sha(HERE/'ranks-m5.json'),
              'incumbent_count_reproduced':str(incumbent),'incumbent_sha256':sha(options.incumbent),
              'denominator_A2_upper':str(A2_upper),'denominator_sha256':sha(options.denominator),
              'rank_gate_permutations':exhaustive,'interleaving_gates_per_quadrant_word_sizes':[1,2,3],
              'all_safe_roots':summary(book.m,safe,A2_upper,incumbent),
              'paired_seam_roots_also_allowed':summary(book.m,larger,A2_upper,incumbent)}
    with options.output.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    for key in ('all_safe_roots','paired_seam_roots_also_allowed'):
        for row in result[key]['thresholds']:
            print(key,row['strict_manhattan_threshold'],row['original_conditional_lower_decimal'],row['factor_over_incumbent'])


if __name__ == '__main__':
    main()
