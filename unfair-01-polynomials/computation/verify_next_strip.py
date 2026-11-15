#!/usr/bin/env python3
"""Independent integer Sturm arithmetic for the m=2r+2 finite strip.

No SymPy. Positive primitive pseudo-remainders preserve Sturm signs.
This verifier reconstructs polynomials from the recurrence, not stored strings.
"""
from fractions import Fraction as Q
from functools import reduce
from math import gcd
import argparse
import hashlib
import platform
import json
import time
from pathlib import Path


def trim(p):
    p = list(p)
    while len(p)>1 and p[-1] == 0: p.pop()
    return p


def add(a,b):
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                 for i in range(max(len(a),len(b)))])


def scale(p,c): return trim([c*x for x in p])


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)


def primitive(p):
    p=trim(p);g=reduce(gcd,p,0)
    return [x//g for x in p] if g else p


def prem(a,b):
    # Each reduction multiplies by a positive number, then subtracts a
    # multiple of b. Thus the final result is a positive scalar times rem(a,b).
    p=primitive(a);b=primitive(b)
    assert b != [0]
    while p != [0] and len(p)>=len(b):
        shift=len(p)-len(b);c=p[-1]*(1 if b[-1]>0 else -1)
        p=scale(p,abs(b[-1]))
        for j,x in enumerate(b):p[j+shift]-=c*x
        p=primitive(p)
    return p


def sturm(p):
    p=primitive(p);assert len(p)>1
    seq=[p,primitive([i*p[i] for i in range(1,len(p))])]
    while seq[-1] != [0]:
        r=scale(prem(seq[-2],seq[-1]),-1)
        if r==[0]:break
        seq.append(r)
    return seq


def sign_at(p,x):
    x=Q(x);a,b=x.numerator,x.denominator
    value=p[-1];den=1
    for c in reversed(p[:-1]):
        den*=b;value=value*a+c*den
    return (value>0)-(value<0)


def variations(seq,x):
    signs=[v for p in seq if (v:=sign_at(p,x))]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def count(seq,lo,hi):
    assert lo<hi and sign_at(seq[0],lo) and sign_at(seq[0],hi)
    return variations(seq,lo)-variations(seq,hi)


def U(n):return [(-1)**i for i in range(n+1)]

def W(n):return [(i+1)*(-1)**i for i in range(n+1)]

D=[1,-2,3];DD=mul(D,D)


def numerators(r):
    T=U(r+2)
    F=add(add(mul(U(2*r+1),DD),scale(mul(mul(T,W(r+1)),D),-1)),mul(mul(T,T),[1,-3]))
    return {'u1':F,'u0':add(F,[0]+DD)}


def witness(r,which):
    T=U(r+2)
    if which=='overflow':
        return add(add(mul([0,-1,2],DD),mul(mul(T,add(U(r+4),scale([0]+U(r+3),-1))),D)),
                   mul(mul(T,T),add([0]+W(3),scale(W(4),-1))))
    k=int(which)
    assert r<=k<2*r
    return add(mul(D,U(k)),scale(mul(T,W(k-r)),-1))


def self_test():
    for p,lo,hi,expected in [([-1,0,1],-2,2,2),([1,-2,1],0,2,1),
                             ([0,-1,0,1],-2,2,3),([1,0,1],-2,2,0)]:
        assert count(sturm(p),Q(lo),Q(hi))==expected
        assert count(sturm(scale(p,-7)),Q(lo),Q(hi))==expected


def strip_endpoint_zeros(p,lo,hi):
    p=trim(p);removed={str(lo):0,str(hi):0}
    for point in (lo,hi):
        assert point in (0,1) or sign_at(p,point)
        while not sign_at(p,point):
            if point==0:
                assert p[0]==0;p=p[1:]
            else:
                quotient=[0]*(len(p)-1);quotient[-1]=p[-1]
                for k in range(len(quotient)-2,-1,-1):quotient[k]=p[k+1]+quotient[k+1]
                assert p[0]+quotient[0]==0;p=quotient
            removed[str(point)]+=1
    return p,removed


def main():
    ap=argparse.ArgumentParser();ap.add_argument('ledger',type=Path)
    ap.add_argument('--small-ledger',type=Path,default=Path(__file__).resolve().parent/'results/011-small-sturm-ledger.json')
    ap.add_argument('--last',type=int,default=71)
    args=ap.parse_args();assert 9<=args.last<=71
    self_test();rows=json.loads(args.ledger.read_text());byr={row['r']:row for row in rows}
    assert len(byr)==len(rows)
    small=json.loads(args.small_ledger.read_text());assert [row['r'] for row in small]==list(range(3,9))
    for row in small:
        byr[row['r']]={'r':row['r'],'branches':{
            key:{'sturm_roots_open_interval':entry['strict_open_root_count'],
                 'isolators':[{'lo':root['I'][0],'hi':root['I'][1]} for root in entry['roots']]}
            for key,entry in row['branches'].items()}}
    results=[];started=time.monotonic()
    for r in range(3,args.last+1):
        lo,hi=(Q(0),Q(1)) if r<9 else (Q(13,20),Q(17,20))
        row=byr[r];branches={}
        for name,p in numerators(r).items():
            p,removed=strip_endpoint_zeros(p,lo,hi)
            seq=sturm(p);number=count(seq,lo,hi)
            entry=row['branches'][name]
            assert number==entry['sturm_roots_open_interval']
            roots=[];previous=lo
            for index,interval in enumerate(entry.get('isolators',[])):
                a,b=Q(interval['lo']),Q(interval['hi'])
                assert previous<a<b<hi and count(seq,a,b)==1
                previous=b
                chosen='overflow' if r==9 else str(r+4 if r==11 and index==1 else r)
                wp=witness(r,chosen);ws=sturm(wp)
                assert count(ws,a,b)==0
                sign=sign_at(wp,(a+b)/2)
                assert sign==(1 if chosen=='overflow' else -1)
                roots.append({'interval':[str(a),str(b)],'witness':chosen,'sign':sign})
            assert len(roots)==number
            branches[name]={'count':number,'variations':[variations(seq,lo),variations(seq,hi)],
                            'excluded_endpoint_multiplicities':removed,'roots':roots}
        results.append({'r':r,'search_interval':[str(lo),str(hi)],'branches':branches})
        print(json.dumps({'verified_r':r,'elapsed_seconds':time.monotonic()-started}),flush=True)
    output={'range':[3,args.last],'all_3_to_71_verified':args.last==71,'results':results,
            'python':platform.python_version(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'ledger_sha256':hashlib.sha256(args.ledger.read_bytes()).hexdigest(),
            'small_ledger_sha256':hashlib.sha256(args.small_ledger.read_bytes()).hexdigest()}
    out=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round011')
    assert Path('/Volumes/data').is_mount();out.mkdir(parents=True,exist_ok=True)
    (out/'independent-sturm.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
