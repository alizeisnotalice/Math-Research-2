#!/usr/bin/env python3
"""Exact active-basis recovery for nine36-m16 line 8 from the owned integer cache."""
from __future__ import annotations
import hashlib,json,subprocess,sys
from itertools import combinations_with_replacement,combinations,product
from pathlib import Path
import numpy as np
import sympy as S
from scipy.linalg import qr
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent; CACHE=HERE/'line8-d2a2-matrixcache'; TOPIC=Path('/Volumes/data/CodexProjects/T2/unfair-01-polynomials')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pair(x):x=S.Rational(x);return [int(x.p),int(x.q)]
def main():
 assert Path('/Volumes/data').is_mount()
 capture,candidate,receipt=[HERE/x for x in ('line8-active-basis-capture.json','line8-recovered-candidate.json','line8-independent-verification.json')]
 assert not any(x.exists() for x in (capture,candidate,receipt))
 m=json.loads((CACHE/'matrix-metadata.json').read_text()); A=np.load(CACHE/m['cache_file'],mmap_mode='r')
 assert sha(CACHE/m['cache_file'])==m['cache_file_sha256'] and hashlib.sha256(np.asarray(A).tobytes()).hexdigest()==m['integer_matrix_sha256']
 nleft,natom=m['free_multiplier_columns'],m['positive_atom_columns']; b=np.zeros(A.shape[0]);b[-1]=1
 ans=linprog(np.r_[np.ones(2*nleft),np.zeros(natom)],A_eq=A,b_eq=b,bounds=(0,None),method='highs');assert ans.success,ans.message
 active=np.flatnonzero(ans.x>1e-9);raw={'format':'round034-line8-active-basis-v1','cache_metadata_sha256':sha(CACHE/'matrix-metadata.json'),'matrix_sha256':m['integer_matrix_sha256'],'solver':{'message':ans.message,'objective':float(ans.fun)},'active_columns':[{'index':int(k),'float_hex':float(ans.x[k]).hex()}for k in active]}
 capture.write_text(json.dumps(raw,ensure_ascii=False,indent=2)+'\n')
 _,r,piv=qr(np.asarray(A[:,active],dtype=float).T,pivoting=True,mode='economic');rank=int(np.sum(abs(np.diag(r))>1e-9));rows=[int(k)for k in piv[:rank]];raw['recovery']={'active_count':len(active),'numeric_rank':rank,'selected_rows':rows}
 if rank!=len(active):raw['recovery']['status']='dependent_active_support';capture.write_text(json.dumps(raw,ensure_ascii=False,indent=2)+'\n');print(json.dumps(raw['recovery']));return
 exact=[S.Rational(x)for x in S.Matrix(np.asarray(A[rows][:,active],dtype=np.int64).tolist()).inv()*S.Matrix([int(b[x])for x in rows])]
 good=(S.Matrix(np.asarray(A[:,active],dtype=np.int64).tolist())*S.Matrix(exact)==S.Matrix([int(x)for x in b]) and all(x>=0 for x in exact));raw['recovery'].update({'exact_full_matrix_nonnegative':good,'min_value':str(min(exact))})
 if not good:raw['recovery']['status']='exact_active_basis_failed';capture.write_text(json.dumps(raw,ensure_ascii=False,indent=2)+'\n');print(json.dumps(raw['recovery']));return
 names=m['strict_variables']; monos=[()]+[(x,)for x in names]+list(combinations_with_replacement(names,2)); products=[(k,z)for k in m['equation_diagonals']for z in monos];atoms=[]
 for size in (1,2):
  for picks in combinations(range(len(names)),size):
   for signs in product((0,1),repeat=size):atoms.append((tuple(names[t]for t,s in zip(picks,signs)if s==0),tuple(names[t]for t,s in zip(picks,signs)if s==1)))
 assert len(products)==nleft and len(atoms)==natom
 z=[S.Integer(0)]*A.shape[1]
 for k,v in zip(active,exact):z[int(k)]=v
 p=[z[k]-z[nleft+k]for k in range(nleft)];q=z[2*nleft:]
 rec={'format':'positive-certificate-v1','case':{'name':'nine36-m16-line8-cache-active-basis','I':m['I'],'J':m['J']},'strict_variables':names,'unit_variables':m['unit_variables'],'basis':{'multiplier_total_degree_at_most':2,'positive_atom_total_degree_at_most':2,'atom_form':'product(positive) * product(1-minus)'},'certificate':{'multipliers':[{'diagonal':k,'monomial':list(z),'coefficient':pair(v)}for v,(k,z)in zip(p,products)if v],'positive_atoms':[{'positive':list(pos),'one_minus':list(neg),'coefficient':pair(v)}for v,(pos,neg)in zip(q,atoms)if v],'exact_recovery_check':True},'status':'certificate_found'}
 candidate.write_text(json.dumps(rec,ensure_ascii=False,indent=2)+'\n');run=subprocess.run([sys.executable,str(TOPIC/'computation/verify_positive_certificate.py'),str(candidate)],text=True,capture_output=True);receipt.write_text(json.dumps({'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr},ensure_ascii=False,indent=2)+'\n');raw['recovery'].update({'independent_fraction_verifier_returncode':run.returncode,'status':'certificate_verified'if run.returncode==0 else'candidate_rejected'});capture.write_text(json.dumps(raw,ensure_ascii=False,indent=2)+'\n');print(json.dumps(raw['recovery'],ensure_ascii=False))
if __name__=='__main__':main()
