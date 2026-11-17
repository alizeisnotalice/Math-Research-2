#!/usr/bin/env python3
"""Exact finite audit of the partial-coverage interface on all 208 subsets."""
import argparse, hashlib, json, sys, time
from pathlib import Path
P,FULL=64,(1<<64)-1
def parity(x):return x.bit_count()&1
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def rref(rows,d=15):
 p=[0]*d
 for x in rows:
  for b in range(d-1,-1,-1):
   if x>>b&1:
    if p[b]:x^=p[b]
    else:
     for c in range(b-1,-1,-1):
      if x>>c&1 and p[c]:x^=p[c]
     p[b]=x
     for c in range(b+1,d):
      if p[c]>>b&1:p[c]^=x
     break
 return p
def nullspace(rows,d=15):
 p=rref(rows,d);out=[]
 for f in range(d):
  if not p[f]:
   x=1<<f
   for b in range(d):
    if p[b] and parity(p[b]&x):x|=1<<b
   out.append(x)
 return out
def span(gens):
 a=[0]
 for g in gens:a += [x^g for x in a]
 return a
def combine(gens,z):
 out=0
 for j,g in enumerate(gens):
  if z>>j&1:out^=g
 return out
def fwht(a):
 a=a[:];s=1
 while s<P:
  for i in range(0,P,2*s):
   for j in range(i,i+s):a[j],a[j+s]=a[j]+a[j+s],a[j]-a[j+s]
  s*=2
 return a
def qmask(perm,h):return sum((((((x&7)&perm[x>>3]).bit_count()^((h>>(x>>3))&1))&1)<<x) for x in range(P))
def dirs(mask):
 pts=[x for x in range(P) if mask>>x&1];b=pts[0];return {x^b for x in pts}
def basis(S):
 p={}
 for x in S:
  y=x
  while y:
   b=y.bit_length()-1
   if b in p:y^=p[b]
   else:p[b]=y;break
 return list(p.values())
def affine_rows(mask,ds,G):
 b=(mask&-mask).bit_length()-1;out=[]
 for u in ds:
  row=0
  for j,g in enumerate(G):
   sl=((g>>(b^u))^(g>>b))&1
   if any((((g>>(x^u))^(g>>x))&1)!=sl for x in range(P) if mask>>x&1):raise AssertionError('nonaffine')
   row|=sl<<j
  out.append(row)
 return out
def components(G):
 words=span(G);mw=[x for x in words if x.bit_count()==8];adj={i:{j for j,y in enumerate(mw) if i!=j and (mw[i]&y).bit_count()==4} for i in range(len(mw))};seen=set();out=[]
 for i in range(len(mw)):
  if i in seen:continue
  todo=[i];seen.add(i);c=[]
  while todo:
   u=todo.pop();c.append(u)
   for v in adj[u]:
    if v not in seen:seen.add(v);todo.append(v)
  if len(c)==6:
   T=0
   for j in c:T|=mw[j]
   U=set.intersection(*(dirs(mw[j]) for j in c))
   if T.bit_count()!=16 or len(U)!=4 or rank(U)!=2:raise AssertionError('E4')
   out.append(('E4',T,basis(U)))
  elif len(c)==1:
   B=mw[c[0]];D=dirs(B)
   if len(D)!=8 or rank(D)!=3:raise AssertionError('iso')
   out.append(('iso',B,basis(D)))
  else:raise AssertionError('component')
 return out,words
ap=argparse.ArgumentParser();ap.add_argument('--candidates',type=Path,required=True);ap.add_argument('--extensions',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();beg=time.monotonic()
rc,re,ra=a.candidates.read_bytes(),a.extensions.read_bytes(),a.audit.read_bytes();CAND,EXT,AUD=json.loads(rc),json.loads(re),json.loads(ra)
seed=next(r for r in AUD['rows'] if r['source_file'].endswith('corrected-results.json') and r['row_index']==1);q=qmask(seed['permutation'],seed['h_bits']);W=CAND['W_generators'];allrows=[];enum_count=0;hier=[]
for ri,rec in enumerate(EXT['records']):
 G=W+rec['quotient_basis'];comps,words=components(G);h=sum(t=='E4' for t,_,_ in comps);p=sum(t=='iso' for t,_,_ in comps)
 if rank(G)!=15 or (h,p) not in ((4,0),(2,4)):raise AssertionError('G')
 C=[FULL]+[sum(((x>>j)&1)<<x for x in range(P)) for j in range(6)]
 hier.append({'record_index':ri,'dim_V_mod_C_minwords':15-rank(C+[x for x in words if x.bit_count()==8]),'dim_V_mod_C_weight_le_16':15-rank(C+[x for x in words if x.bit_count()<=16])})
 for sb in range(1<<len(comps)):
  chosen=[comps[j] for j in range(len(comps)) if sb>>j&1];m=sum(2 if t=='E4' else 1 for t,_,_ in chosen);used=0;L=[]
  for t,mask,ds in chosen:
   if used&mask:raise AssertionError('overlap')
   used|=mask;L+=affine_rows(mask,ds,G)
  R=FULL^used;hh=sum(t=='E4' for t,_,_ in chosen);pp=sum(t=='iso' for t,_,_ in chosen);RL=rank(L);kb=nullspace(L)
  if len(kb)!=15-RL or any(parity(row&z) for row in L for z in kb):raise AssertionError('kernel')
  kw=[combine(G,z) for z in span(kb)];rwords=[x&R for x in kw]
  gp={x.bit_length()-1:x for x in basis(G)}
  for value in set(rwords):
   reduced=value
   while reduced:
    bit=reduced.bit_length()-1
    if bit not in gp:raise AssertionError('residual zero-extension not in V')
    reduced^=gp[bit]
  dwr=rank(rwords);covered_kernel=[x for x in kw if not (x&R)];covered_rank=rank([x&used for x in covered_kernel]);budget=RL+3*hh+pp+dwr
  if 15!=budget or len(kb)-dwr!=3*hh+pp or covered_rank!=3*hh+pp or dwr>(8-m)**2:raise AssertionError(('exact budget',ri,sb,RL,len(kb),dwr,hh,pp,covered_rank,budget))
  e={'record_index':ri,'subset_bits':sb,'h':hh,'p':pp,'m':m,'R_size':R.bit_count(),'rankL':RL,'dimK':len(kb),'dimV_R':dwr,'covered_kernel_dimension':len(kb)-dwr,'covered_pattern_rank':covered_rank,'dimension_budget_rhs':budget,'all_residual_zero_extensions_in_V':True,'R_minimum_word_count':sum((x&~R)==0 for x in words if x.bit_count()==8)}
  if not R:
   if dwr:raise AssertionError('emptyR')
   e.update({'R_empty':True,'atoms':0,'all_K_residual_support_enumerated':False});allrows.append(e);continue
  DR={i for i,x in enumerate(fwht([(-1 if q>>y&1 else 1) if R>>y&1 else 0 for y in range(P)])) if x}
  if len(DR)!=R.bit_count():raise AssertionError('DR')
  indep=[]
  for x in rwords:
   if rank(indep+[x])>len(indep):indep.append(x)
  atoms={}
  for x in range(P):
   if R>>x&1:
    lab=tuple((v>>x)&1 for v in indep);atoms[lab]=atoms.get(lab,0)|(1<<x)
  if dwr>len(atoms) or any(x.bit_count()*R.bit_count()<64 for x in atoms.values()):raise AssertionError('atoms')
  for atom in atoms.values():
   supp={i for i,x in enumerate(fwht([(-1 if q>>y&1 else 1) if atom>>y&1 else 0 for y in range(P)])) if x}
   if not supp<=DR:raise AssertionError('atom support')
  exhaustive=sb.bit_count()==1
  if exhaustive:
   for v in kw:
    output=fwht([(-1 if (q^v)>>y&1 else 1) if R>>y&1 else 0 for y in range(P)])
    if any(x not in (-8,0,8) for x in output):raise AssertionError('residual amplitude')
    supp={i for i,x in enumerate(output) if x}
    if supp!=DR:raise AssertionError('K support')
   enum_count+=len(kw)
  e.update({'R_empty':False,'seed_residual_Walsh_support_size':len(DR),'atoms':len(atoms),'min_atom_times_R':min(x.bit_count()*R.bit_count() for x in atoms.values()),'all_atom_supports_in_fixed_D_R':True,'all_K_residual_support_enumerated':exhaustive,'K_slices_enumerated':len(kw) if exhaustive else 0});allrows.append(e)
if len(allrows)!=208:raise AssertionError(len(allrows))
out={'status':'passed','scope':'finite all-subset interface audit; full K enumeration only marked singleton subsets','rows':allrows,'weight_hierarchy':hier,'all_K_slices_enumerated_total':enum_count,'elapsed_seconds':time.monotonic()-beg,'candidate_sha256':hashlib.sha256(rc).hexdigest(),'extensions_sha256':hashlib.sha256(re).hexdigest(),'audit_sha256':hashlib.sha256(ra).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version}
a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','subsets':len(allrows),'all_K_slices':enum_count,'elapsed_seconds':out['elapsed_seconds']}))
