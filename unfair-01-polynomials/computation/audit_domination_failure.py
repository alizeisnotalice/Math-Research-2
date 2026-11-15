#!/usr/bin/env python3
"""Exhaust all injections and issue Hall obstructions independently of the matcher."""
import itertools
import json
from collections import defaultdict


def audit(I,J):
    edges=defaultdict(list)
    for i,j in itertools.product(I,J): edges[i+j].append((i,j))
    ua={i for i in I if i in (0,max(I)) or len(edges[i])==1}
    ub={j for j in J if j in (0,max(J)) or len(edges[j])==1}
    def monomial(edge):
        i,j=edge
        return ({('a',i)} if i not in ua else set())|({('b',j)} if j not in ub else set())
    terms={k:[monomial(e) for e in es] for k,es in edges.items()}
    certificates=[];injections=0;weak_pairs=0;hall_pairs=0
    for source,target in itertools.permutations(sorted(edges),2):
        P,Q=terms[source],terms[target]
        weak=[]
        for match in itertools.permutations(range(len(P)),len(Q)):
            injections+=1
            if all(P[match[t]]<=Q[t] for t in range(len(Q))):
                assert len(P)==len(Q) and all(P[match[t]]==Q[t] for t in range(len(Q))), 'Strict domination found'
                weak.append(list(match))
        if weak:
            weak_pairs+=1
            certificates.append({'source':source,'target':target,'weak_only_matchings':weak})
        else:
            hall_pairs+=1
            witness=None
            for size in range(1,len(Q)+1):
                for subset in itertools.combinations(range(len(Q)),size):
                    neighbors={p for t in subset for p in range(len(P)) if P[p]<=Q[t]}
                    if len(neighbors)<len(subset):
                        witness={'target_subset':list(subset),'source_neighbors':sorted(neighbors)}
                        break
                if witness: break
            assert witness is not None
            certificates.append({'source':source,'target':target,'hall':witness})
    assert ua=={0,max(I)} and ub=={0,max(J)}
    return {'I':I,'J':J,'multiplicities':{k:len(edges[k]) for k in sorted(edges)},
            'unit_I':sorted(ua),'unit_J':sorted(ub),'pairs':len(certificates),
            'injections_tested':injections,'hall_pairs':hall_pairs,'weak_only_pairs':weak_pairs,
            'strict_domination_exists':False,'certificates':certificates}


if __name__=='__main__':
    I=[0,3,5,7];J=[0,3,5,6,8,9,11,15,17]
    result=[audit(I,J),audit(I,sorted(J+[13]))]
    print(json.dumps({'cases':result},indent=2))
