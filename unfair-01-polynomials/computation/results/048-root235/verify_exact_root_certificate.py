#!/usr/bin/env python3
"""Independent exact-Fraction checker for make_exact_root_certificate output.

It does not import the generator.  Every Rouché inequality, root-disk
separation, conjugacy/group constraint, frozen mask mapping, interval factor,
constant normalization, and recorded negative upper endpoint is recomputed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(0)


def U(x): return F(int(x['n']), int(x['d']))
def P(x): return {'n': str(x.numerator), 'd': str(x.denominator)}
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def l1(x): return abs(x[0])+abs(x[1])
def li(x): return max(abs(x[0]),abs(x[1]))


def poly(c,z):
    q=(F(0),F(0))
    for a in reversed(c): q=add(mul(q,z),(F(a),F(0)))
    return q


def deriv(c,z):
    q=(F(0),F(0))
    for k in range(len(c)-1,0,-1): q=add(mul(q,z),(F(k*c[k]),F(0)))
    return q


def remainder(c,z,r):
    rho=l1(z); total=F(0)
    for k,a in enumerate(c):
        if k>=2 and a: total += abs(a)*((rho+r)**k-rho**k-k*r*rho**(k-1))
    return total


class I:
    __slots__=('m','e')
    def __init__(self,m,e=F(0)):
        self.m,self.e=F(m),F(e)
        if self.e<0: raise AssertionError('negative interval radius')
    def plus(self,y): return I(self.m+y.m,self.e+y.e)
    def times(self,y): return I(self.m*y.m,abs(self.m)*y.e+abs(y.m)*self.e+self.e*y.e)
    def inv(self):
        lo,hi=self.m-self.e,self.m+self.e
        if lo<=0<=hi: raise AssertionError('normalization interval crosses zero')
        u,v=F(1,lo),F(1,hi); lo,hi=min(u,v),max(u,v)
        return I((lo+hi)/2,(hi-lo)/2)
    def upper(self): return self.m+self.e


def pmul(a,b):
    out=[I(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]=out[i+j].plus(x.times(y))
    return out


def factor(group,r):
    if len(group)==1:
        a,b=group[0]; assert b==0
        return [I(-a,r),I(1)]
    assert len(group)==2 and group[1]==(group[0][0],-group[0][1])
    a,b=group[0]; return [I(a*a+b*b,2*(abs(a)+abs(b))*r+r*r),I(-2*a,2*r),I(1)]


def bit_vector(g,mask):
    assert 0<=mask<(1<<(g-1))
    return [0]+[(mask>>(g-1-i))&1 for i in range(1,g)]


def normal(factors,select):
    raw=[I(1)]
    for b,f in zip(select,factors):
        if b: raw=pmul(raw,f)
    inv=raw[0].inv()
    return [I(1)]+[v.times(inv) for v in raw[1:]],raw[0]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--certificate',type=Path,required=True)
    ap.add_argument('--input',type=Path,required=True); ap.add_argument('--diagnostics',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True); args=ap.parse_args()
    cert=json.loads(args.certificate.read_text()); inp=json.loads(args.input.read_text()); diag=json.loads(args.diagnostics.read_text())
    assert cert['input_sha256']==hashlib.sha256(args.input.read_bytes()).hexdigest()
    assert cert['diagnostics_sha256']==hashlib.sha256(args.diagnostics.read_bytes()).hexdigest()
    r=U(cert['radius']); assert r>0
    declared_ids=[int(i) for i in cert['tasks']]
    actual_ids=[int(row['task']) for row in cert['certificates']]
    assert declared_ids==sorted(set(declared_ids))
    assert actual_ids==declared_ids and len(actual_ids)==len(set(actual_ids))
    assert declared_ids==[18,19,20,21,22,23]
    wanted={int(x['id']):x for x in inp['tasks'] if int(x['id']) in declared_ids}
    assert sorted(wanted)==declared_ids
    expected={i:[] for i in wanted}
    for row in diag:
        if row.get('task') in expected and row.get('classification')=='nearunknown': expected[row['task']].append(int(row['mask']))
    checked=[]; total=0; least_margin=None; least_upper=None
    for tc in cert['certificates']:
        tid=int(tc['task']); task=wanted[tid]; c=[int(v) for v in task['C']]
        assert all(v in (0,1) for v in c) and len(c)-1==tc['degree']
        assert task['groups']==tc['groups']==len(tc['centres_by_group'])
        groups=[]
        for raw,source_group in zip(tc['centres_by_group'],task['roots']):
            assert len(raw)==len(source_group)
            group=[(U(z['re']),U(z['im'])) for z in raw]
            if len(group)==1: assert group[0][1]==0
            elif len(group)==2: assert group[1]==(group[0][0],-group[0][1])
            else: raise AssertionError('bad group cardinality')
            # Bind every recorded exact disk centre to its original binary64 seed;
            # this does not prove the seed, but prevents silently changing group order.
            for (a,b), seed in zip(group,source_group):
                sr, si=F(str(seed[0])),F(str(seed[1]))
                assert abs(a-sr)<F(1,10**8) and abs(b-si)<F(1,10**8)
            groups.append(group)
        flat=[z for g in groups for z in g]; assert len(flat)==len(c)-1
        assert min(max(abs(a[0]-b[0]),abs(a[1]-b[1])) for i,a in enumerate(flat) for b in flat[:i])>2*r
        for z in flat:
            left=li(deriv(c,z))*r; right=l1(poly(c,z))+remainder(c,z,r)
            assert left>right
            margin=left-right; least_margin=margin if least_margin is None or margin<least_margin else least_margin
        factors=[factor(g,r) for g in groups]
        rows=tc['mask_certificates']; assert [int(x['mask']) for x in rows]==expected[tid]
        for row in rows:
            mask=int(row['mask']); bits=bit_vector(task['groups'],mask)
            assert bits==row['bits_by_group']
            a,a0=normal(factors,bits); b,b0=normal(factors,[1-x for x in bits])
            assert a0.m-a0.e>0 or a0.m+a0.e<0
            assert b0.m-b0.e>0 or b0.m+b0.e<0
            assert row['negative_factor'] in ('A','B')
            seq=a if row['negative_factor']=='A' else b
            idx=int(row['negative_index']); assert 0<=idx<len(seq)
            upper=seq[idx].upper(); assert upper<0 and P(upper)==row['negative_upper']
            least_upper=upper if least_upper is None or upper>least_upper else least_upper
            total+=1
        checked.append({'task':tid,'masks':len(rows),'roots':len(flat)})
    assert total==sum(len(v) for v in expected.values())
    report={'status':'passed_independent_exact_fraction_recheck','certificate':str(args.certificate),
            'certificate_sha256':hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
            'tasks':checked,'masks':total,'radius':P(r),'least_rouche_margin':P(least_margin),
            'least_negative_upper':P(least_upper),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'masks':total,'out':str(args.out)}))


if __name__=='__main__': main()
