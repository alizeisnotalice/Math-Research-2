#!/usr/bin/env python3
"""Bounded exact/interval tests of R57 algebra; no geometric enumeration."""
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from math import isqrt
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode=True
path=Path(__file__).with_name('r057_quadratic_endpoint_kernel.py');spec=importlib.util.spec_from_file_location('r057kernel',path);k=importlib.util.module_from_spec(spec);spec.loader.exec_module(k)
begin=time.process_time();checks={}

def sqrtbox(q,bits):
    assert q>=0;D=1<<bits;n=isqrt((q.numerator<<(2*bits))//q.denominator);lo=F(n,D)
    return lo,lo if lo*lo==q else F(n+1,D)

def kbox(x,bits):
    lo,hi=sqrtbox(x.field.Q,bits)
    return (x.a+x.b*lo,x.a+x.b*hi) if x.b>=0 else (x.a+x.b*hi,x.a+x.b*lo)

def oracle_field_sign(x):
    if not x.a and not x.b:return 0
    for bits in (16,32,64,128,256,512,1024):
        lo,hi=kbox(x,bits)
        if lo>0:return 1
        if hi<0:return -1
    raise AssertionError('field oracle interval unresolved')

def ebox(e,bits):
    a,b=kbox(e.c,bits)
    if not e.s:return a,b
    u,v=kbox(e.U,bits)
    if u<0:return None
    lo=sqrtbox(u,bits)[0];hi=sqrtbox(v,bits)[1]
    return (a+lo,b+hi) if e.s>0 else (a-hi,b-lo)

def oracle_endpoint_sign(a,b):
    for bits in (32,64,128,256,512,1024):
        x,y=ebox(a,bits),ebox(b,bits)
        if x is None or y is None:continue
        if x[0]>y[1]:return 1
        if x[1]<y[0]:return -1
    raise AssertionError('endpoint oracle unresolved; do not infer equality from interval overlap')

f=k.Field(2);E=f.element
# Known conjugate squares: sqrt(3+2sqrt2)=1+sqrt2;
# sqrt(3-2sqrt2)=sqrt2-1, whose positive branch is essential.
edge=[(k.Endpoint(E(0),1,E(3,2)),k.Endpoint(E(2),1,E(3,-2)),0),
      (k.Endpoint(E(-1),1,E(2)),k.Endpoint(E(0),1,E(3,-2)),0),
      (k.Endpoint(E(1),-1,E(2)),k.Endpoint(E(0),1,E(3,-2)),-1),
      (k.Endpoint(E(0),1,E(3,2)),k.Endpoint(E(2),-1,E(3,-2)),1),
      (k.Endpoint(E(0),1,E(0)),k.Endpoint(E(0),-1,E(0)),0),
      (k.Endpoint(E(-3),1,E(4)),k.Endpoint(E(0),0,E(0)),-1),
      (k.Endpoint(E(-1),1,E(4)),k.Endpoint(E(0),0,E(0)),1)]
for a,b,expected in edge:
    assert k.compare_endpoints(a,b)==expected and k.compare_endpoints(b,a)==-expected
checks['branch_checked_edge_cases']=len(edge)*2
for Q in [F(0),F(4),F(1,9),F(2),F(3,5),F(1)+F(1,2**200)]:
    field=k.Field(Q)
    for a in [-2,-1,0,1,2]:
        for b in [-2,-1,0,1,2]:
            x=field.element(a,b);assert x.sign()==oracle_field_sign(x)
            if x.sign():assert x/x==1 and (1/x)*x==1
checks['field_sign_and_inverse_grid']=150
# Independent exact oracle whenever U=p^2: the principal sqrt is |p| in K.
rootcases=0
for Q in [F(2),F(3,5)]:
    field=k.Field(Q);rows=[]
    for aa in [-2,-1,0,1,2]:
      for bb in [-2,-1,0,1,2]:
        p=field.element(aa,bb);absolute=p if oracle_field_sign(p)>=0 else -p
        for c in [F(-1),F(0),F(1)]:
          for s in [-1,1]:rows.append((k.Endpoint(field.element(c),s,p*p),field.element(c)+s*absolute))
    for i,(left,knownleft) in enumerate(rows):
        for shift in (1,37,91):
            right,knownright=rows[(i+shift)%len(rows)]
            expected=oracle_field_sign(knownleft-knownright)
            assert k.compare_endpoints(left,right)==expected;rootcases+=1
checks['quadratic_square_radical_comparisons']=rootcases
# Generic nested roots: independently separate rational interval enclosures.
event=E(F(2,5),F(-1,7));generic=[]
for c in [F(-2),F(-1,3),F(0),F(1,5),F(2)]:
 for q in [F(1),F(2),F(5)]:
  for slope in [F(-1,3),F(0),F(1,4)]:
   for s in [-1,1]:generic.append(k.recipe_at_event([c,0,s,q,slope],0,event))
genericcases=0
for i,left in enumerate(generic):
 for shift in (1,17):
    right=generic[(i+shift)%len(generic)];expected=oracle_endpoint_sign(left,right)
    assert k.compare_endpoints(left,right)==expected;genericcases+=1
checks['generic_nested_radical_interval_comparisons']=genericcases
# Distinguish true crossing, second-order tangency, mass-persistent but
# different full two-parameter recipes, and square-q0 nonpersistent ties.
zero=E(0)
p=lambda row:k.prepare_recipe(row,0,zero)
cross=k.contact(p([0,0,1,1,1]),p([1,0,0,0,0]));assert cross['contact_order']==1 and cross['left_sign']==-1 and cross['right_sign']==1
tangent=k.contact(p([0,0,1,1,1]),p([-1,0,1,4,2]));assert tangent['contact_order']==2 and tangent['left_sign']==tangent['right_sign']==-1
persistent=k.contact(p([0,0,1,4,0]),p([2,1,0,0,0]));assert persistent['persistent_in_t'] and not persistent['same_full_recipe']
assert k.persistent_recipe_key([0,0,1,4,1])!=k.persistent_recipe_key([2,0,0,0,0])
checks['contact_order_controls']={'crossing':cross,'tangency':tangent,'t_persistent_z_distinct':persistent}
# Affine-radicand origin changes must preserve identities and event values.
r1=[0,0,1,2,3];r2=[0,0,1,F(7,2),3]
assert k.persistent_recipe_key(r1,0)==k.persistent_recipe_key(r2,F(1,2))
assert k.compare_endpoints(k.recipe_at_event(r1,0,event),k.recipe_at_event(r2,F(1,2),event))==0
checks['parameter_origin_identity']=True
rejected=[]
for name,fn in [('negative_Q',lambda:k.Field(-1)),('negative_radicand',lambda:k.Endpoint(E(0),1,E(-1))),
 ('zero_division',lambda:E(1)/E(0)),('nonanalytic_zero_mass_jet',lambda:k.jet_from_recipe([0,0,1,0,1],0,zero,1))]:
 try:fn()
 except (ValueError,ZeroDivisionError):rejected.append(name)
 else:raise AssertionError('unsafe case accepted '+name)
checks['fail_closed_cases']=rejected
result={'status':'passed','scope':'bounded exact quadratic/nested-radical and first/second jet algebra tests; no R57 geometry enumeration',
 'checks':checks,'kernel_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-begin,'random_seed':None,'floating_tolerance':None,'actual_command':[sys.executable,*sys.argv]}
if not os.path.ismount('/Volumes/data'):raise RuntimeError('/Volumes/data is not mounted')
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
