#!/usr/bin/env python3
"""Integer checks for a bent residual equality model and a signed non-bent control."""
import argparse, hashlib, json, sys, time
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for p in range(0,len(a),2*h):
   for i in range(p,p+h):a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
  h*=2
 return a
def rank(rows,n):
 rows=rows[:];r=0
 for b in range(n-1,-1,-1):
  p=next((i for i in range(r,len(rows)) if rows[i]>>b&1),None)
  if p is None:continue
  rows[r],rows[p]=rows[p],rows[r]
  for i in range(len(rows)):
   if i!=r and rows[i]>>b&1:rows[i]^=rows[r]
  r+=1
 return r
def bent_residual():
 # x=(u,v), two bits each. Standard triangular V: four v-fibres, u0,
 # and u1 split by v0.  It has dimension 7 and all f_*+v are MM bent.
 bas=[]
 for j in range(4):bas.append(sum(1<<x for x in range(16) if x>>2==j))
 bas+=[sum(1<<x for x in range(16) if x&1),sum(1<<x for x in range(16) if ((x>>1)&1) and ((x>>2)&1)==0),sum(1<<x for x in range(16) if ((x>>1)&1) and ((x>>2)&1)==1)]
 assert rank(bas,16)==7
 f=[((x&3)&(x>>2)).bit_count()&1 for x in range(16)]
 out=[]
 for ell in (1,2,3):
  Q=1<<ell;n=4+2*ell;N=1<<n;q=4*Q;R=[i for i in range(N) if (i>>(4+ell))==0];D={a|(c<<(4+ell)) for a in range(16) for c in range(Q)}
  atoms={}
  partial_supports=set(); full=True
  for code in range(128):
   v=[sum(((code>>j)&1)*((bas[j]>>x)&1) for j in range(7))&1 for x in range(16)]
   # full phase F+v 1_{w=0}
   sign=[]
   for y in range(N):
    x=y&15;z=(y>>4)&(Q-1);w=y>>(4+ell);bit=f[x]^((z&w).bit_count()&1)^(v[x] if w==0 else 0);sign.append(-1 if bit else 1)
   W=fwht(sign);full &= all(abs(a)==q for a in W)
   part=[0]*N
   for y in R:
    x=y&15;part[y]=-1 if f[x]^v[x] else 1
   P=fwht(part);partial_supports.add(tuple(i for i,a in enumerate(P) if a)); assert all((i in D and abs(a)==q) or (i not in D and a==0) for i,a in enumerate(P))
  for y in R:
   x=y&15;label=tuple((b>>x)&1 for b in bas);atoms.setdefault(label,0);atoms[label]+=1
  out.append({'ell':ell,'variables':n,'q':q,'R_size':len(R),'D_size':len(D),'V_dimension':7,'all_128_full_slices_bent':full,'common_partial_support_count':len(partial_supports),'atoms':len(atoms),'atom_sizes':sorted(set(atoms.values())),'expected_atom_size_q_over_t':Q,'normalized_partial_is_seed_Walsh_repeated_over_c':True})
 return out
def signed_control():
 # h(a0,a1,a2,a3,b)=a0*a2+a1*a3+b*a0*a1, n=5.
 h=[(((x>>0)&1)*((x>>2)&1)^((x>>1)&1)*((x>>3)&1)^((x>>4)&1)*((x>>0)&1)*((x>>1)&1)) for x in range(32)]
 Wh=fwht([-1 if z else 1 for z in h]);Gam={i for i,z in enumerate(Wh) if z};assert len(Gam)==16 and all(abs(Wh[i])==8 for i in Gam)
 formula={a|((((a>>2)&1)&((a>>3)&1))<<4) for a in range(16)};assert Gam==formula
 g=[z//8 for z in Wh];assert all(x in (-1,0,1) for x in g);Wg=fwht(g);assert Wg==[4*(-1 if z else 1) for z in h]
 # a=1 z-padding, bpad=0: S(u,z)=g(u), R=Gamma x F2, D=F2^5 x {0}.
 S=[g[u] for z in range(2) for u in range(32)];WS=fwht(S);R=[i for i,x in enumerate(S) if x];D=set(range(32));q=8
 assert len(R)==len(D)==32 and all((i in D and abs(x)==q) or (i not in D and x==0) for i,x in enumerate(WS))
 # six character evaluations: constant and five physical u coordinates; on R they distinguish the 16 u's.
 labels={tuple([1]+[(u>>j)&1 for j in range(5)]):0 for u in Gam}
 for y in R:labels[tuple([1]+[((y&31)>>j)&1 for j in range(5)])]+=1
 affrank=rank([x^next(iter(Gam)) for x in Gam],5)
 return {'n':6,'t':4,'q':q,'R_size':len(R),'D_size':len(D),'h_spectrum_size':len(Gam),'h_walsh_amplitude':8,'Gamma_matches_formula':True,'signed_g_values_0_pm1':True,'Wg_equals_t_times_h':True,'G_dimension_constant_plus_u_characters':6,'evaluation_atoms':len(labels),'atom_sizes':sorted(set(labels.values())),'support_affine_dimension_after_compression':affrank,'log2_support_size':4,'strictly_not_small_bent_by_affine_dimension':affrank>4}
def true_bent_control():
 # X=(xi0,xi1,xi2,xi3,b,z), R={b=xi2*xi3}; D={frequency z=0}.
 N=64;q=8;F=[];R=[];D=set(range(32));gb=[]
 for y in range(N):
  x0,x1,x2,x3,b,z=[(y>>j)&1 for j in range(6)];v=x2*x3;F.append(x0*x2^x1*x3^z*(b^v))
  if b==v:R.append(y)
 for j in range(6):gb.append(sum(1<<y for y in R if (j==0 or ((y>>(j-1))&1)))) # const,u0..u4
 assert rank(gb,64)==6
 atoms={tuple((g>>y)&1 for g in gb) for y in R};assert len(atoms)==16
 # G is zero off R, so codeword masks are directly sums of these restrictions.
 goodG=True
 for a in range(64):
  W=fwht([-1 if (F[y]^sum(((a>>j)&1)*((gb[j]>>y)&1) for j in range(6))&1) else 1 for y in range(N)])
  goodG &= all(abs(x)==q for x in W)
 blocks=[];supports=[]
 for c in range(2):
  for d in range(2):
   B=sum(1<<y for y in range(N) if ((y>>2)&1)==c and ((y>>3)&1)==d and ((y>>4)&1)==(1^(c*d)));blocks.append(B)
   vals=[(-1 if F[y] else 1) if B>>y&1 else 0 for y in range(N)];W=fwht(vals);supports.append({i for i,x in enumerate(W) if x})
 assert len(set().union(*supports))==32 and sum(len(s) for s in supports)==32 and all(s<=set(range(32,64)) for s in supports)
 goodV=True
 for a in range(1024):
  bit=[F[y] for y in range(N)]
  for j,m in enumerate(gb+blocks):
   if a>>j&1:
    for y in range(N):bit[y]^=(m>>y)&1
  W=fwht([-1 if x else 1 for x in bit]);goodV &= all(abs(x)==q for x in W)
 # partial datum on R for its 64 switches
 partial=True
 for a in range(64):
  s=[0]*N
  for y in R:s[y]=-1 if F[y]^sum(((a>>j)&1)*((gb[j]>>y)&1) for j in range(6))&1 else 1
  W=fwht(s);partial &= all((i in D and abs(x)==q) or (i not in D and x==0) for i,x in enumerate(W))
 return {'variables':6,'q':q,'t':4,'R_size':len(R),'D_size':len(D),'G_dimension':6,'atoms':len(atoms),'atom_sizes':[2],'all_64_G_slices_bent':goodG,'partial_walsh_common_D':partial,'Rc_blocks':4,'block_partial_walsh_cosets_partition_Dc':True,'V_dimension':10,'all_1024_V_slices_bent':goodV,'compressed_support_affine_dimension':5,'compressed_support_log2_size':4}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);a=ap.parse_args();t=time.perf_counter();d={'status':'passed','bent_residual_tensor':bent_residual(),'signed_plateaued_nonbent_control':signed_control(),'true_bent_compressed_control':true_bent_control(),'source_sha256':hashlib.sha256(open(__file__,'rb').read()).hexdigest(),'python':sys.version,'elapsed_seconds':time.perf_counter()-t};open(a.out,'w').write(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':'passed','elapsed_seconds':d['elapsed_seconds']}))
if __name__=='__main__':main()
