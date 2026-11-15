#!/usr/bin/env python3
"""Bounded numerical adversary; no negative search claim is a proof."""
import ctypes,json,math,os,sys,time,hashlib,platform
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
LIB='/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r053_cell_union.dylib'
L=ctypes.CDLL(LIB);D=ctypes.c_double*9;I=ctypes.c_int*9
for name in ('full','heuristic','area'):getattr(L,name).restype=ctypes.c_double
OUT=Path(__file__).resolve().parent/'results/r053-mass-adversary'

def rows(w):return [sum(w[:3]),sum(w[6:]),sum(w[0::3]),sum(w[2::3])]
def value(w,seeds=5,full=False):
 p=I();a=L.full(D(*w),p,0) if full else L.heuristic(D(*w),p,seeds,530015)
 return a,list(p)
def main():
 assert os.path.ismount('/Volumes/data');OUT.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter();rng=np.random.default_rng(530015);records=[];counts={'sample':0,'feasible':0,'optimization':0}
 def record(w,tag,seeds=5):
  a,p=value(w,seeds);r={'w':list(map(float,w)),'area':a,'order':p,'boundary_rows':rows(w),'tag':tag};records.append(r);return a
 for alpha in (.08,.2,.5,1,2,5):
  for _ in range(500):
   w=rng.dirichlet([alpha]*9);counts['sample']+=1
   if min(rows(w))>.25:counts['feasible']+=1;record(w,'dirichlet-'+str(alpha))
 for w in ([.25,0,.25,0,0,0,.25,0,.25],[.125,.125,.125,.125,0,.125,.125,.125,.125],[1/9]*9):record(w,'baseline',20)
 sampled=sorted(records,key=lambda r:-r['area'])[:12]
 def obj(x):
  w=np.asarray(x)/sum(x);side=rows(w);counts['optimization']+=1
  if min(side)<.25:return 100+100*(.25-min(side))
  return -value(w,5)[0]
 opt=[]
 for seed in (530015,530016,530017):
  res=differential_evolution(obj,[(1e-6,1)]*9,seed=seed,popsize=8,maxiter=90,tol=1e-7,polish=False,workers=1,updating='immediate')
  w=res.x/sum(res.x);record(w,'DE-'+str(seed),50);opt.append({'seed':seed,'success':bool(res.success),'message':str(res.message),'nit':res.nit,'nfev':res.nfev})
 leaders=sorted(records,key=lambda r:-r['area'])[:15]
 for r in leaders:r['complete_permutation_min'],r['complete_best_order']=value(r['w'],full=True)
 result={'scope':'Floating bounded maximization of best-order nine-cell union area; only numerical evidence. Full enumeration of 9! orders only at reported leader points. Empty cells skipped.','seed':530015,'counts':counts,'sample_leaders':sampled,'optimization':opt,'leaders':leaders,'elapsed_seconds':time.perf_counter()-start,'python':sys.version,'numpy':np.__version__,'scipy':__import__('scipy').__version__,'platform':platform.platform(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'kernel_sha256':hashlib.sha256(Path(__file__).with_name('r053_cell_union.cpp').read_bytes()).hexdigest()}
 (OUT/'bounded-search.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
