#!/usr/bin/env python3
"""Search and check exact termwise domination witnesses on fixed positive supports."""
import argparse
import json
from pathlib import Path


def find_certificate(I,J):
    edges={}
    for i in I:
        for j in J:
            edges.setdefault(i+j,[]).append((i,j))
    labels=[('a',i) for i in I]+[('b',j) for j in J]
    strict=[label for label in labels if label[1] not in (0,(I if label[0]=='a' else J)[-1]) and len(edges[label[1]])>1]
    bit={label:1<<k for k,label in enumerate(strict)}
    masks={k:[bit.get(('a',i),0)|bit.get(('b',j),0) for i,j in es] for k,es in edges.items()}
    for source in sorted(edges):
        P=masks[source]
        for target in sorted(edges):
            Q=masks[target]
            if source==target or len(P)<len(Q): continue
            candidates=[[i for i,p in enumerate(P) if p&q==p] for q in Q]
            matched={}
            def augment(t,seen):
                for si in candidates[t]:
                    if si in seen: continue
                    seen.add(si)
                    if si not in matched or augment(matched[si],seen):
                        matched[si]=t
                        return True
                return False
            if not all(augment(t,set()) for t in sorted(range(len(Q)),key=lambda x:len(candidates[x]))): continue
            mapping={t:si for si,t in matched.items()}
            if len(P)>len(Q) or any(P[mapping[t]]!=Q[t] for t in range(len(Q))):
                return {'source':source,'target':target,'target_to_source':[mapping[t] for t in range(len(Q))]}
    return None


def verify_certificate(I,J,certificate):
    # Recompute each requested diagonal by membership, not the searcher's masks.
    def diagonal(k): return [(i,k-i) for i in I if k-i in set(J)]
    def factors(edge):
        i,j=edge
        factors=set()
        for name,index,support in [('a',i,I),('b',j,J)]:
            multiplicity=len(diagonal(index))
            if index not in (0,support[-1]) and multiplicity>1:
                factors.add((name,index))
        return factors
    P=diagonal(certificate['source']); Q=diagonal(certificate['target'])
    mapping=certificate['target_to_source']
    assert P and Q and len(mapping)==len(Q) and len(set(mapping))==len(mapping)
    assert all(isinstance(i,int) and 0<=i<len(P) for i in mapping)
    proper=False
    for edge,si in zip(Q,mapping):
        p,q=factors(P[si]),factors(edge)
        assert p<=q
        proper |= p<q
    assert len(P)>len(Q) or proper, 'No strict inequality certified'
    return True


def self_check():
    cases=[([0,1],[0,2],False),([0,1],[0,1],True),
           ([0,1,2,4,9],[0,1,3,4,6,11],True),
           ([0,1,2,4,7],[0,1,2,4,8,10,13],True),
           ([0,1,2,3,5,11],[0,1,3,4,5,7,13],True),
           ([0,1,2,4,6,11],[0,1,3,4,6,8,13],True)]
    for I,J,expected in cases:
        cert=find_certificate(I,J)
        assert (cert is not None)==expected,(I,J,cert)
        if cert: verify_certificate(I,J,cert)
    try:
        verify_certificate([0,1],[0,2],{'source':0,'target':3,'target_to_source':[0]})
    except AssertionError:
        pass
    else:
        raise AssertionError('Accepted non-strict equality as contradiction')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('unknown',type=Path)
    args=ap.parse_args()
    self_check()
    output=[]
    for line in args.unknown.read_text().splitlines():
        c=json.loads(line);cert=find_certificate(c['I'],c['J'])
        if cert: verify_certificate(c['I'],c['J'],cert)
        output.append({**c,'certificate':cert,'status':'excluded' if cert else 'unknown'})
    print(json.dumps({'self_checks_passed':True,'cases':output},indent=2))
