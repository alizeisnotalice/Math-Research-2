"""Exact all-word marked-path and reverse-slot audit."""
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import permutations
from math import comb,perm
from pathlib import Path
import argparse
import json
import sys

D=Path(__file__).resolve().parent


def table(path):
    return {tuple(int(v) for v in row.split()[:4]):int(row.split()[4])
            for row in path.read_text().splitlines()}


def brute(base):
    sys.path.insert(0,str(base))
    from audit_distance_layers import grid,rho,roots
    from verify_contact import trace
    n=3
    adj=grid(n)
    totals=Counter()
    counts=Counter()
    variants=defaultdict(Counter)
    lengths={}
    for word in permutations(range(9)):
        peaks=roots(word,adj)
        if len(peaks)!=2:
            continue
        p,q=peaks
        rho0=rho(p,q,n)
        totals[rho0]+=1
        if rho0>=3:
            continue
        reduced,j0,K,U=trace(word,q,adj)
        if j0>6:
            continue
        valid={v for v in U if rho(p,v,n)>=3 or rho(q,v,n)%2}
        m=9-j0-1
        slots={v:i for i,v in enumerate(reduced[j0:],1)}
        qslot=word.index(q)
        def walk(at,seen,path):
            if path and rho(p,at,n)>=3:
                yield path
            if len(path)<4:
                for v in valid-seen:
                    if adj[at] & (1<<v):
                        yield from walk(v,seen|{v},path+(v,))
        for path in walk(q,set(),()):
            k=len(path)
            y=path[-1]
            used=set(path)
            new_red=K+tuple(reversed((q,)+path[:-1]))+tuple(v for v in reduced[j0:] if v not in used)
            out=new_red[:qslot]+(y,)+new_red[qslot:]
            assert roots(out,adj)==(p,y)
            assert trace(out,y,adj)[1:3]==(j0,K)
            prefix=new_red[j0:j0+k]
            mark=lambda v: rho(p,v,n)<3 and rho(q,v,n)%2==0
            assert next(v for v in new_red[j0:] if mark(v))==q
            assert prefix[-1]==q
            parity=(q%n+q//n)%2
            key=(out,parity)
            assert lengths.setdefault(key,k)==k
            counts[j0+1,rho0,3,k]+=1
            variants['all'][key]+=1
            selected=[slots[v] for v in path]
            max_slot=min(m,k+1)
            if selected==sorted(selected) and selected[-1]<=max_slot:
                variants['ordered_early'][key]+=1
                if selected[0]==1:
                    variants['ordered_first'][key]+=1
    bounds={}
    for typ,cases in variants.items():
        load=Fraction(0)
        max_load=Fraction(0)
        for key,mult in cases.items():
            out,parity=key
            p,y=roots(out,adj)
            m=9-trace(out,y,adj)[1]-1
            k=lengths[key]
            maximum=min(m,k+1)
            cap={'all':perm(m,k),'ordered_early':comb(maximum,k),
                 'ordered_first':comb(maximum-1,k-1)}[typ]
            assert mult<=cap and rho(p,y,n)>=3
            pressure=Fraction(mult,cap)
            load+=pressure
            max_load=max(max_load,pressure)
        scale=sum(totals.values())
        bounds[typ]={'images':len(cases),'averaged_probability_lower':str(load/(2*scale)),
                     'max_individual_image_weight':str(max_load),
                     'total_path_word_pairs':sum(cases.values())}
    return counts,bounds


def main():
    p=argparse.ArgumentParser()
    p.add_argument('base',type=Path)
    args=p.parse_args()
    observed,b=brute(args.base)
    assert observed==table(D/'n3.tsv')
    exact=json.loads((args.base/'contact-verification.json').read_text())
    rows=[]
    for c in exact['sizes']:
        n=c['n']
        T=c['A2']
        paths=table(D/f'n{n}.tsv')
        for threshold in sorted({t for j,d,t,k in paths}):
            weighted=sum(Fraction(count,2*perm(n*n-j,k)) for (j,d,t,k),count in paths.items() if t==threshold)
            near=next(z['source_below_target'] for z in c['reachability'] if z['target_at_least']==threshold)
            far=T-near
            assert 0<weighted<=far
            rows.append({'n':n,'target_at_least':threshold,'paths':sum(x for key,x in paths.items() if key[2]==threshold),
                         'certified_probability_lower':str(weighted/T),'decimal_probability_lower':float(weighted/T),
                         'actual_tail':str(Fraction(far,T))})
    print(json.dumps({'n3_exact_capacity_audit':b,'all_path_bounds':rows},indent=2))


if __name__=='__main__':
    main()
