#!/usr/bin/env python3
"""Exact FWHT verification of the R31 scalar rank-four gamma-patch example."""
from __future__ import annotations
import hashlib,json,platform,sys,time
from pathlib import Path
OUT=Path(__file__).with_name("r31_scalar_rank4_absorption_result.json")
N=8;QX=1<<N;q=16
def par(x):return x.bit_count()&1
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for s in range(0,len(a),2*h):
   for j in range(h):
    left,right=a[s+j],a[s+j+h];a[s+j],a[s+j+h]=left+right,left-right
  h*=2
 return a

def rank(vs):
 b={}
 for v in vs:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return len(b)
def indep(vs):
 out=[];r=0
 for v in vs:
  if rank(out+[v])>r:out.append(v);r+=1
 return out
def F(x):
 x0,x1,x2,x3,x4,x5,u,v=((x>>i)&1 for i in range(8))
 P=x1*x2;Q=P^(x3*x4);return P^(x0*(x1^x3))^((x5^Q)*(x2^x4))^(u*v)
def Qfun(x):return ((x>>1)&1)*((x>>2)&1)^(((x>>3)&1)*((x>>4)&1))
def block(a,b):return sum(1<<x for x in range(QX) if ((x>>0)&1)==a and ((x>>5)&1)==(Qfun(x)^b) and not((x>>7)&1))
def affine(j,x):return 1 if j==0 else (x>>(j-1))&1
def lift(fun,x):return fun(x) if not((x>>7)&1) else 0
def word(fun):return sum(lift(fun,x)<<x for x in range(QX))
def span(gens):
 s={0}
 for g in gens:s|={x^g for x in list(s)}
 return s
def normal_basis(direction):
 return indep([n for n in range(QX) if all(par(n&d)==0 for d in direction)])
def affine_direction(supp):
 base=supp[0];return base,tuple(sorted(x^base for x in supp))
def local_quad(fun,a,b):
 # ANF quadratic bits on y=(x1,x2,x3,x4,u); Mobius transform of the five-variable restriction.
 table=[]
 for y in range(32):
  x=(a<<0)|(((y>>0)&1)<<1)|(((y>>1)&1)<<2)|(((y>>2)&1)<<3)|(((y>>3)&1)<<4)|((Qfun((y<<1))^b)<<5)|(((y>>4)&1)<<6)
  table.append((fun>>x)&1)
 coef=table[:]
 for bit in range(5):
  for m in range(32):
   if m>>bit&1:coef[m]^=coef[m^(1<<bit)]
 return sum(coef[(1<<i)|(1<<j)]<<(i*5+j) for i in range(5) for j in range(i+1,5))
def main():
 started=time.perf_counter()
 blocks=[block(a,b) for a in (0,1) for b in (0,1)]
 assert sum(x.bit_count() for x in blocks)==128 and not any(blocks[i]&blocks[j] for i in range(4) for j in range(i))
 # Seven old global affine functions plus four block indicators; retain an independent rank-nine basis.
 candidates=[word(lambda x,j=j:affine(j,x)) for j in range(7)]+blocks
 gens=indep(candidates);assert len(gens)==9
 assert all(B.bit_count()==32 for B in blocks)
 R=0
 for B in blocks:R|=B
 assert all(not(g & (((1<<QX)-1)^R)) for g in gens)
 code_to_param={}
 for p in range(512):
  w=0
  for j,g in enumerate(gens):
   if p>>j&1:w^=g
  code_to_param[w]=p
 blockparams=[code_to_param[b] for b in blocks]
 offsets=span(blockparams);assert rank(blockparams)==4
 # The global x5 affine direction is retained as the seventh independent candidate.
 tau=1<<6
 assert gens[6]==candidates[6]
 D=frozenset(xi for xi in range(QX) if not((xi>>6)&1))
 tile_dirs=None;fixed=None;labels=[];allbent=partialok=True
 local_quads=[]
 for a in (0,1):
  for b in (0,1):local_quads.append([local_quad(g,a,b) for g in gens])
 for p in range(512):
  signs=[]
  for x in range(QX):
   bit=F(x)
   for j,g in enumerate(gens):
    if p>>j&1:bit^=((g>>x)&1)
   signs.append(1 if not bit else -1)
  allbent &= all(abs(v)==q for v in fwht(signs))
  coverage=[0]*QX;current_dirs=[];current_labels=[]
  for bi,B in enumerate(blocks):
   part=fwht([signs[x] if B>>x&1 else 0 for x in range(QX)])
   supp=[xi for xi,v in enumerate(part) if v]
   tauval=(p>>6)&1
   for gamma in (0,1):
    tile=[xi for xi in supp if (((xi>>5)&1)^tauval)==gamma]
    partialok &= len(tile)==q and all(abs(part[xi])==q for xi in tile)
    for xi in tile:coverage[xi]+=1
    base,direction=affine_direction(tile);current_dirs.append(direction)
    normals=normal_basis(direction);assert len(normals)==4
    current_labels.extend(par(n&base) for n in normals)
  union=frozenset(i for i,c in enumerate(coverage) if c)
  partialok &= all(c==1 for c in coverage if c) and union==D
  if fixed is None:fixed,tile_dirs=union,current_dirs
  else:partialok &= fixed==union and tile_dirs==current_dirs
  labels.append(sum(bit<<i for i,bit in enumerate(current_labels)))
 rel=[x^labels[0] for x in labels];kernel={p for p,x in enumerate(rel) if not x}
 for p in range(512):
  expected=0
  for j in range(9):
   if p>>j&1:expected^=rel[1<<j]
  assert rel[p]==expected
 # Every block sees only the scalar Q quadratic direction; modulo <Q> the freeze rank is zero.
 qrank_checks=[]
 for qs in local_quads:
  expected_Q=(1<<1)|(1<<13)
  qrank_checks.append(qs==[expected_Q if j==6 else 0 for j in range(9)])
 result={"scope":"Exact verification of the displayed scalar rank-four absorption example only.","python":sys.version,"platform":platform.platform(),"q":q,"base_rank4_graph_blocks":4,"code_rank":len(gens),"phase_count":512,"all_phases_bent":bool(allbent),"residual_size":128,"M":8,"block_size":32,"gamma_tiles":8,"gamma_tile_size":q,"fixed_D_size":len(fixed),"fixed_D_is_u_frequency_zero":fixed==D,"partial_gamma_tiles_fixed_disjoint_pm_q":bool(partialok),"joint_labelled_translation_rank":rank(rel),"translation_labels_linear_all_parameters":True,"translation_kernel_equals_block_constants":kernel==offsets,"block_constant_kernel_rank":rank(blockparams),"per_block_quadratic_variation_scalar_Q":all(qrank_checks),"freeze_rank_mod_scalar_Q":0,"tau_x5_parameter_mask":tau,"elapsed_seconds":time.perf_counter()-started,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");print(json.dumps(result,ensure_ascii=False,indent=2))
 assert allbent and partialok and rank(rel)==5 and kernel==offsets and all(qrank_checks)
if __name__=="__main__":main()
