#!/usr/bin/env python3
"""Independent numerical root-pair/layout verifier for parameterized coprime caches."""
from __future__ import annotations
import argparse, hashlib, itertools, json
from math import comb
from pathlib import Path
import mpmath as mp
import numpy as np
mp.mp.dps=100
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def conv(a,b):
 z=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]+=x*y
 return z
def stretch(h,q):z=[0]*((len(h)-1)*q+1);z[::q]=h;return z
def franel(n,q):return sum(comb(q,k)**n for k in range(q+1))
def local(a,b,u,v):
 z=[]
 for ca in itertools.combinations(range(a),u):
  for cb in itertools.combinations(range(b),v):z.append((sum(1<<x for x in tuple(ca)+tuple(a+j for j in cb)),tuple(ca),tuple(cb)))
 return sorted(z)
def unpair(row):return np.asarray(row['hi'],dtype=np.float64)+np.asarray(row['lo'],dtype=np.float64)
def cmplx(x):return mp.mpc(x[0],x[1])
def resid(h,z):return abs(sum(c*z**i for i,c in enumerate(h)))/sum(abs(c)*abs(z)**i for i,c in enumerate(h))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--options',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);x=ap.parse_args()
 raw=json.loads(x.input.read_text())['candidates'][0]; cache=json.loads(x.options.read_text()); h=list(raw['H_coefficients_ascending']);a,b=int(raw['a']),int(raw['b']);n=(len(h)-1)//2;terminal=a*b
 assert len(h)%2 and h[0]==h[-1]==1 and raw['kind']==cache['kind']=='coprime_substitution'
 assert (cache['candidate'],cache['a'],cache['b'],cache['n_blocks'],cache['balance_terminal_E'])==(raw['name'],a,b,n,terminal)
 C=conv(stretch(h,a),stretch(h,b)); ch=hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest();assert C==raw['C_coefficients_ascending'] and ch==raw['C_sha256'] and set(C)<={0,1}
 assert len(raw['blocks'])==len(cache['blocks'])==raw['parent_block_count']==n and len(raw['groups'])==raw['total_groups']==n*(a+b)
 g={int(q['index']):q for q in raw['groups']};assert len(g)==len(raw['groups']); max_parent=mp.mpf(0);max_lift=mp.mpf(0);max_quad=0.0
 for rb in raw['blocks']:
  alpha=cmplx(rb['parent_root_positive_imag']);pr=resid(h,alpha);max_parent=max(max_parent,pr);assert mp.im(alpha)>0 and pr<mp.mpf('1e-75')
  for family,power,ids in [('a_lift',a,rb['a_lift_groups']),('b_lift',b,rb['b_lift_groups'])]:
   assert len(ids)==power
   for j,i in enumerate(ids):
    q=g[int(i)];assert (q['block'],q['family'],q['lift_index'])==(rb['block'],family,j) and q['parent_root_positive_imag']==rb['parent_root_positive_imag']
    z=cmplx(q['lift_root']);lr=abs(z**power-alpha)/max(1,abs(alpha));max_lift=max(max_lift,lr);assert lr<mp.mpf('1e-75')
    expect=[float(abs(z)**2),float(-2*mp.re(z)),1.0];err=max(abs(r-s) for r,s in zip(q['coefficients_ascending'],expect));max_quad=max(max_quad,err);assert err<5e-15
 balance={}
 for u in range(a+1):
  for v in range(b+1):
   E=b*u+a*v
   if 0<E<terminal:assert E not in balance;balance[E]=(u,v)
 es=sorted(balance);assert [z['E'] for z in raw['E_options']]==es+[terminal]; layouts={}
 for rb,cb in zip(raw['blocks'],cache['blocks']):
  assert rb['block']==cb['block'];ai=list(rb['a_lift_groups']);bi=list(rb['b_lift_groups']);assert cb['a_lift_groups']==ai and cb['b_lift_groups']==bi and rb['all_group_indices']==ai+bi==cb['all_group_indices'];assert sorted(map(int,cb['options']))==es+[terminal]
  for E in es:
   u,v=balance[E];rows=cb['options'][str(E)];want=local(a,b,u,v);assert len(rows)==len(want)==comb(a,u)*comb(b,v);layouts.setdefault(E,[]).append(len(rows))
   for r,(mask,ca,bb) in zip(rows,want):
    selected=[ai[j] for j in ca]+[bi[j] for j in bb];assert r['mode']=='mixed' and (r['u'],r['v'])==(u,v) and r['selected_local_mask_decimal']==str(mask) and r['selected_group_indices']==selected and [e['group_index'] for e in r['selected_local_pairs']]==selected
    assert len(r['Ablock']['hi'])==len(r['Ablock']['lo'])==2*(u+v)+1 and len(r['Bblock']['hi'])==len(r['Bblock']['lo'])==2*(a+b-u-v)+1
  tr=cb['options'][str(terminal)];assert [r['mode'] for r in tr]==['all_a','all_b'] and tr[0]['selected_group_indices']==ai and tr[1]['selected_group_indices']==bi;layouts[terminal]=[2]*n
 internal=sum(int(np.prod(layouts[E])) for E in es);expected_internal=franel(n,a)*franel(n,b)//2-2;assert internal==expected_internal==raw['balanced_enumeration']['internal_total'];total=internal+2**(n-1);expected=franel(n,a)*franel(n,b)//2+2**(n-1)-2;assert total==expected==raw['balanced_enumeration']['nonconstant_complement_representatives']
 A=np.array([1.]);B=np.array([1.])
 for cb in cache['blocks']:
  fair=cb['options'][str(terminal)][0];A=np.convolve(A,unpair(fair['Ablock']));B=np.convolve(B,unpair(fair['Bblock']))
 A/=A[0];B/=B[0]
 err_a=float(np.max(np.abs(A-np.asarray(stretch(h,a),dtype=np.float64))))
 err_b=float(np.max(np.abs(B-np.asarray(stretch(h,b),dtype=np.float64))))
 ferr=max(err_a,err_b);assert ferr<1e-9
 report={'status':'passed_numerical_rootpair_input_and_layout_checks_only','candidate':raw['name'],'a':a,'b':b,'n_blocks':n,'integer_binary_C':True,'block_count':n,'total_groups':len(raw['groups']),'E_options':es+[terminal],'radices_by_E':{str(E):layouts[E] if E!=terminal else [1]+[2]*(n-1) for E in es+[terminal]},'canonical_total':total,'franel_formula_total':expected,'max_parent_root_relative_residual':mp.nstr(max_parent,10),'max_lift_power_relative_residual':mp.nstr(max_lift,10),'max_quadratic_coeff_float64_error':max_quad,'cached_fair_all_a_max_abs_error':ferr,'input_sha256':sha(x.input),'options_sha256':sha(x.options),'verifier_sha256':sha(Path(__file__)),'claim_boundary':'Numerical root/cache validation only; no exact sign exclusion or counterexample is certified.'}
 x.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
