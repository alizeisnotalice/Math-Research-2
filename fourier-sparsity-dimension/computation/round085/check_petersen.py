"""Exact graphic, cut, signed-Walsh and cyclic-packing certificate."""
from fractions import Fraction
from itertools import combinations
import json
from circuit_fields import basis


def main():
    edges = [(i,(i+1)%5) for i in range(5)] + [(5+i,5+(i+2)%5) for i in range(5)] + [(i,i+5) for i in range(5)]
    cols = [(1<<i) | (1<<j) for i,j in edges]
    assert len(basis(cols)) == 9
    cuts = [sum(((t>>i)^(t>>j)) & 1 for i,j in edges) for t in range(1,1023)]
    assert min(cuts)==3
    eligible = []
    for k in (10,11):
        hist = [0]*15
        n = 0
        for ix in combinations(range(15),k):
            selected = [cols[i] for i in ix]
            if len(basis(selected))!=9 or any(len(basis(selected[:j]+selected[j+1:]))<9 for j in range(k)):
                continue
            n += 1
            for i in ix:
                hist[i] += 1
        assert max(hist)==min(hist)
        eligible.append(dict(edges=k,subgraphs=n,incidence=hist[0]))
    assert eligible[0]['subgraphs']==0 and eligible[1]['subgraphs']>0
    chosen = [cols[i] for i in (0,1,2,3,4,5,6,7,8,11,14)]
    assert all(len(basis(chosen[:j]+chosen[j+1:]))==9 for j in range(11))
    table = [(-1)**(sum((x >> (2*i) & 3)==3 for i in range(7))%2) for x in range(1<<14)]
    for b in range(14):
        step = 1<<b
        for j in range(0,len(table),2*step):
            for i in range(j,j+step):
                a,c = table[i],table[i+step]
                table[i],table[i+step] = a+c,a-c
    assert set(map(abs,table))=={128}
    assert Fraction(20,16*2**14) <= Fraction(11,15*16*128)
    return dict(labels=len(edges), minimum_cut=min(cuts), proper_vertex_sets=len(cuts), eligible=eligible,
                dense_walsh_terms=len(table), rank=28, sparsity=16*(15+2**14), L='2/3', cyclic='11/15')


if __name__=='__main__':
    print(json.dumps(main(),indent=2))
