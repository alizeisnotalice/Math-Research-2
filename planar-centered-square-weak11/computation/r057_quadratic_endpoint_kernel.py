#!/usr/bin/env python3
"""Exact R57 quadratic-event endpoint comparison and one-sided jet ordering.

No floating point or all-pairs comparison. The base ordered field is Q(sqrt(Q));
endpoints are c+s*sqrt(U), c,U in that field and U>=0. Squaring is used only
AFTER opposing branch signs are checked. Recipe rows are [a,v,s,q0,q1], with
radicand q0+q1*(t-t0); a already includes the fixed physical z shift. v is kept
as provenance for full two-parameter identity, not differentiated in t.
"""
import argparse
from fractions import Fraction as F
from functools import cmp_to_key
import gzip
import hashlib
import json
from math import isqrt
import os
from pathlib import Path
import platform
import sys
import time


def sgn(x):return (x>0)-(x<0)

class Field:
    __slots__=('Q','rational_root')
    def __init__(self,Q):
        self.Q=F(Q)
        if self.Q<0:raise ValueError('real quadratic field requires Q>=0')
        n,d=isqrt(self.Q.numerator),isqrt(self.Q.denominator)
        self.rational_root=F(n,d) if n*n==self.Q.numerator and d*d==self.Q.denominator else None
    def element(self,a=0,b=0):
        if isinstance(a,Element):
            if a.field.Q!=self.Q or b:raise ValueError('different base fields or nonzero extra coefficient')
            return a
        a,b=F(a),F(b)
        if self.rational_root is not None:a,b=a+b*self.rational_root,F(0)
        return Element(self,a,b)

class Element:
    __slots__=('field','a','b')
    def __init__(self,field,a,b):self.field,self.a,self.b=field,a,b
    def _q(self,x):return self.field.element(x)
    def __add__(self,x):x=self._q(x);return self.field.element(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return self.field.element(-self.a,-self.b)
    def __sub__(self,x):return self+-self._q(x)
    def __rsub__(self,x):return self._q(x)+-self
    def __mul__(self,x):
        x=self._q(x);return self.field.element(self.a*x.a+self.b*x.b*self.field.Q,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self._q(x);den=x.a*x.a-x.b*x.b*self.field.Q
        if not den:raise ZeroDivisionError('zero element of ordered quadratic field')
        return self*self.field.element(x.a/den,-x.b/den)
    def __rtruediv__(self,x):return self._q(x)/self
    def __pow__(self,n):
        if n==2:return self*self
        if n==0:return self.field.element(1)
        if n==1:return self
        raise ValueError('kernel requires only square, first or zero power')
    def __eq__(self,x):
        try:x=self._q(x)
        except (ValueError,TypeError):return False
        return self.a==x.a and self.b==x.b
    def __hash__(self):return hash((self.field.Q,self.a,self.b))
    def __repr__(self):return 'Element('+str(self.a)+','+str(self.b)+';Q='+str(self.field.Q)+')'
    def sign(self):
        sa,sb=sgn(self.a),sgn(self.b)
        if not sb:return sa
        if not sa:return sb
        if sa==sb:return sa
        return sa*sgn(self.a*self.a-self.b*self.b*self.field.Q)
    def pair(self):return [str(self.a),str(self.b)]


def sign_linear_radical(h,k,U):
    """Exact sign(h+k*sqrt(U)); h,k,U are in the same ordered base field."""
    h=h.field.element(h);k=h.field.element(k);U=h.field.element(U)
    su,sh,sk=U.sign(),h.sign(),k.sign()
    if su<0:raise ValueError('negative real radicand')
    if not su or not sk:return sh
    if not sh:return sk
    if sh==sk:return sh
    # Only opposing signs reach squaring; equality here is a true equality.
    return sh*(h*h-k*k*U).sign()

class Endpoint:
    __slots__=('c','s','U','identity')
    def __init__(self,c,sign,U,identity=None):
        if not isinstance(c,Element):raise TypeError('offset must be a base-field element')
        self.c,self.s,self.U,self.identity=c,int(sign),c.field.element(U),identity
        if self.s not in (-1,0,1):raise ValueError('root sign must be -1,0,1')
        if self.s and self.U.sign()<0:raise ValueError('negative endpoint radicand')
        if not self.s:self.U=c.field.element(0)
    def data(self):return {'c':self.c.pair(),'sign':self.s,'U':self.U.pair()}


def compare_endpoints(left,right):
    """Return exact sign(left-right), including genuine event equality."""
    d=left.c-right.c;U,V=left.U,right.U;s,r=left.s,-right.s
    y=sign_linear_radical(d,d.field.element(s),U)
    if not r or not V.sign():return y
    if not y:return r
    if y==r:return y
    # y=d+s sqrt(U) opposes r sqrt(V); sign of total is
    # sign(y)*sign(y^2-V). The latter has only ONE remaining radical.
    return y*sign_linear_radical(d*d+s*s*U-V,2*s*d,U)


def persistent_recipe_key(row,t0=0,include_z_velocity=True):
    a,v,s,q0,q1=row;a,v,q0,q1=map(F,(a,v,q0,q1));s=int(s)
    if s not in (-1,0,1):raise ValueError('invalid recipe root sign')
    if not s:q0=q1=F(0)
    elif not q1:
        if q0<0:raise ValueError('negative constant recipe radicand')
        n,d=isqrt(q0.numerator),isqrt(q0.denominator)
        if n*n==q0.numerator and d*d==q0.denominator:a+=s*F(n,d);s=0;q0=F(0)
    key=(a,v if include_z_velocity else F(0),s,q0-q1*F(t0),q1)
    return key


def recipe_at_event(row,t0,event):
    a,v,s,q0,q1=row;a,q0,q1=map(F,(a,q0,q1));s=int(s)
    U=event.field.element(q0)+q1*(event-F(t0))
    return Endpoint(event.field.element(a),s,U,persistent_recipe_key(row,t0))


def jet_from_recipe(row,t0,event,order):
    if order not in (1,2):raise ValueError('only first and second t-jets are needed')
    e=recipe_at_event(row,t0,event);q1=F(row[4]);zero=event.field.element(0)
    if not e.s or not q1:return Endpoint(zero,0,zero)
    if e.U.sign()<=0:raise ValueError('positive radicand required for analytic jets')
    C=e.s*q1/(2*e.U) if order==1 else -e.s*q1*q1/(4*e.U*e.U)
    sc=C.sign()
    return Endpoint(zero,sc,C*C*e.U)


def prepare_recipe(row,t0,event):
    return (recipe_at_event(row,t0,event),jet_from_recipe(row,t0,event,1),jet_from_recipe(row,t0,event,2),
            persistent_recipe_key(row,t0),persistent_recipe_key(row,t0,False))


def compare_prepared(left,right,side=0):
    center=compare_endpoints(left[0],right[0])
    if center or not side:return center
    for order in (1,2):
        delta=compare_endpoints(left[order],right[order])
        if delta:return delta*(-1 if side<0 and order==1 else 1)
    if left[4]!=right[4]:raise ArithmeticError('unexpected nonpersistent contact of order>=3')
    return 0


def contact(left,right):
    center=compare_endpoints(left[0],right[0])
    if center:return {'center_sign':center,'contact_order':0,'left_sign':center,'right_sign':center}
    for order in (1,2):
        delta=compare_endpoints(left[order],right[order])
        if delta:return {'center_sign':0,'contact_order':order,'first_nonzero_jet_sign':delta,
                         'left_sign':delta*(-1 if order==1 else 1),'right_sign':delta,
                         'persistent_in_t':False,'same_full_recipe':False}
    if left[4]!=right[4]:raise ArithmeticError('unexpected nonpersistent contact of order>=3')
    return {'center_sign':0,'contact_order':None,'left_sign':0,'right_sign':0,
            'persistent_in_t':True,'same_full_recipe':left[3]==right[3]}


def order_cluster(records,event,t0):
    prepared={str(r['id']):prepare_recipe(r['expr'],t0,event) for r in records}
    if len(prepared)!=len(records):raise ValueError('duplicate endpoint IDs')
    ids=list(prepared)
    def groups(side):
        compare=lambda a,b:compare_prepared(prepared[a],prepared[b],side)
        ordered=sorted(ids,key=cmp_to_key(compare));out=[]
        for i in ordered:
            if out and compare(out[-1][-1],i)==0:out[-1].append(i)
            else:out.append([i])
        return out
    center=groups(0);ties=[]
    for group in center:
        if len(group)<2:continue
        # Adjacent comparisons suffice; no all-pairs tie graph is constructed.
        rs=sorted(group,key=cmp_to_key(lambda a,b:compare_prepared(prepared[a],prepared[b],1)))
        ties.append({'members':group,'adjacent_right_order_contacts':[
            {'ids':[a,b],**contact(prepared[a],prepared[b])} for a,b in zip(rs,rs[1:])]})
    return {'center_order_groups':center,'left_order_groups':groups(-1),'right_order_groups':groups(1),'center_ties':ties,'endpoint_count':len(records)}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--cluster-index',type=int);a=p.parse_args()
    begin=time.process_time();raw=a.input.read_bytes();data=json.loads(gzip.decompress(raw) if a.input.suffix=='.gz' else raw)
    field=Field(data['Q']);event=field.element(*data['event']);t0=F(data.get('t0',0));records=data['endpoints'];lookup={str(r['id']):r for r in records}
    clusters=data.get('candidate_clusters',[list(lookup)])
    jobs=range(len(clusters)) if a.cluster_index is None else [a.cluster_index]
    result=[]
    for j in jobs:
        group=[lookup[str(i)] for i in clusters[j]];result.append({'cluster_index':j,**order_cluster(group,event,t0)})
    out={'status':'passed','scope':'exact value/first/second t-jet order inside listed candidate clusters; no global source/bbox completeness or between-cluster separation asserted',
      'Q':str(field.Q),'event':event.pair(),'t0':str(t0),'clusters':result,'input_sha256':hashlib.sha256(raw).hexdigest(),
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),
      'cpu_seconds':time.process_time()-begin,'actual_command':[sys.executable,*sys.argv]}
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','scope','code_sha256','cpu_seconds')}))
if __name__=='__main__':main()
