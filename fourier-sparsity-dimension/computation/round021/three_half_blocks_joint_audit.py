#!/usr/bin/env python3
"""Exact n=4 witness and joint-cylinder audit for R21 three blocks."""
import hashlib,json,platform,sys,time
from pathlib import Path
N=16; FULL=(1<<N)-1
def fwht(a):
 a=list(a);h=1
 while h<N:
  for i in range(0,N,2*h):
   for j in range(i,i+h):x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
  h*=2
 return a
def bent(m):return all(abs(x)==4 for x in fwht([1-2*((m>>z)&1) for z in range(N)]))
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def span(bs):
 S={0}
 for b in bs:S|={x^b for x in tuple(S)}
 return S
def basis(S):
 o=[]
 for x in sorted(S,reverse=True):
  if rank(o+[x])>len(o):o.append(x)
 return o
def affine_flat(S):
 b=next(iter(S));D={x^b for x in S};r=rank(D)
 return len(S)==1<<r and span(basis(D))==D,r
def affine_span_dim(S):
 b=next(iter(S));return rank({x^b for x in S})
def pmask(f):return sum(f(x)<<x for x in range(N))
def partial(f,A):return fwht([(1-2*((f>>x)&1)) if A>>x&1 else 0 for x in range(N)])
def extension(word,A,aff,unique=True):
 hit=[z for z in aff if ((z^word)&A)==0]
 assert len(hit)==1 if unique else len(hit)>=1
 return min(hit)
def enc(v,xi,d):return v|(xi<<d)
def nullspace(rows,n):return [x for x in range(1<<n) if all(((x&r).bit_count()&1)==0 for r in rows)]
def normal_basis(S,n):
 b=next(iter(S));D={x^b for x in S};return basis(nullspace(D,n))
def main():
 t=time.monotonic()
 f0=pmask(lambda x:((x&1)&((x>>1)&1))^(((x>>2)&1)&((x>>3)&1)))
 assert bent(f0)
 legal=[f0^g for g in range(1<<N) if bent(g) and (f0^g).bit_count()==6]
 assert len(legal)==192
 witness=None
 for i,A in enumerate(legal):
  for B in legal[i+1:]:
   if A&B:continue
   if bent(f0^A^B):witness=(A,B);break
  if witness:break
 assert witness
 A,B=witness;R=FULL^A^B
 assert A.bit_count()==B.bit_count()==6 and R.bit_count()==4 and bent(f0^R)
 coords=[pmask(lambda x,j=j:(x>>j)&1) for j in range(4)]
 gens=[FULL]+coords+[A,B];d=rank(gens);assert d==7
 words=[]
 for c in range(1<<d):
  w=0
  for j,g in enumerate(gens):
   if c>>j&1:w^=g
  words.append(w)
 assert len(set(words))==128 and all(bent(f0^w) for w in words)
 aff=[pmask(lambda x,c=c,l=l:c^((x&l).bit_count()&1)) for c in range(2) for l in range(N)]
 affinfo={m:(c,l) for c in range(2) for l in range(N) for m in [pmask(lambda x,c=c,l=l:c^((x&l).bit_count()&1))]}
 blocks=[('A',A),('B',B),('R',R)];cyl=[];direct={};blockinfo=[]
 ambient=d+4
 for name,C in blocks:
  W=partial(f0,C);S={z for z,x in enumerate(W) if x};assert set(W)<={-4,0,4};rr=affine_span_dim(S)
  if C.bit_count()==6:
   assert len(S)==6 and rr==4
   pairs=[{sorted(S)[0],sorted(S)[1]},{sorted(S)[2],sorted(S)[3]},{sorted(S)[4],sorted(S)[5]}]
   dirs=[next(iter(p))^next(iter(p-{next(iter(p))})) for p in pairs]
   assert len(set(dirs))==3
  else:
   ok,_=affine_flat(S);assert len(S)==4 and rr==2 and ok;pairs=[S]
  ext=[];sets=[]
  for vi,w in enumerate(words):
   _,ell=affinfo[extension(w,C,aff,unique=(C.bit_count()==6))];ext.append(ell)
   actual={enc(vi,z,d) for z,x in enumerate(partial(f0^w,C)) if x}
   direct[ name,vi]=actual
  for line in pairs:
   Z={enc(vi,xi,d) for vi,ell in enumerate(ext) for xi in range(N) if (xi^ell) in line}
   ok,rd=affine_flat(Z);assert ok
   sets.append(Z);cyl.append((name,Z,rd,ambient-rd))
  for vi in range(1<<d):
   assert direct[name,vi]==set().union(*[{x for x in Z if (x&((1<<d)-1))==vi} for Z in sets])
  blockinfo.append({'name':name,'mask_hex':hex(C),'weight':C.bit_count(),'baseline_support':sorted(S),'support_affine_dimension':rr,'extension_linear_parts':ext})
 allsets=[x[1] for x in cyl]
 assert all(not(allsets[i]&allsets[j]) for i in range(len(allsets)) for j in range(i))
 norms=[]
 for _,Z,rd,codim in cyl:
  nb=normal_basis(Z,ambient);assert len(nb)==codim;norms+=nb
 pnorm=[x&((1<<d)-1) for x in norms];K=nullspace(pnorm,d)
 kr={}
 for name,C in blocks:
  kr[name]=rank([words[k]&C for k in K]);assert kr[name]<=1
 result={'schema':'a28-r021-three-half-block-joint-v1','scope':'one exact n=4 common-bent witness; finite validation only','arithmetic':'exact F2 and integer FWHT','f0_mask_hex':hex(f0),'legal_weight6_switches':len(legal),'A_hex':hex(A),'B_hex':hex(B),'R_hex':hex(R),'dimension':d,'phase_count':len(words),'all_phases_bent':True,'blocks_individually_legal':True,'blocks':blockinfo,'joint_cylinder_count':len(cyl),'joint_cylinders':[{'block':n,'dimension':rd,'codimension':co} for n,_,rd,co in cyl],'joint_cylinders_disjoint':True,'normal_union_rank':rank(norms),'parameter_projection_rank':rank(pnorm),'K_dimension':rank(K),'K_block_restriction_ranks':kr,'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 out=Path(__file__).with_name('three-half-block-joint-audit.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
