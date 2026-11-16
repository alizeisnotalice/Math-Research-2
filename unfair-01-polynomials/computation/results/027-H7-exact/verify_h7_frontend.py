"""Exact H7(3,4) root coverage and independent real 25-bit-mask ledger.
No coefficient exclusion. Helpers do not import a searcher or its ledger generator.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse,json,hashlib
from h7_exact_common import cj,pf,mother,lift,dis
H=[1,1,1,0,0,0,1,1]
ALL=(1<<25)-1

def require(x,msg):
 if not x:raise ValueError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def local_sets():
 out=[]
 for block in range(3):
  local={}
  for bits in range(128):
   u=(bits&7).bit_count();v=(bits>>3).bit_count();E=4*u+3*v
   local.setdefault(E,[]).append(bits<<(7*block))
  out.append(local)
 real={}
 for bits in range(16):
  e=bits&1;eta=(bits>>1)&1;w=(bits>>2).bit_count();E=4*e+8*eta+6*w
  real.setdefault(E,[]).append(bits<<21)
 out.append(real)
 return out

def layouts():
 sets=local_sets();common=set.intersection(*(set(d) for d in sets));full=[]
 for E in sorted(common):
  full.extend((E,m[0]|m[1]|m[2]|m[3]) for m in product(*(d[E] for d in sets)))
 require(len(full)==len(set(full)),'duplicate full layout')
 fs=set(full);require(all((24-E,ALL^m) in fs for E,m in fs),'complement missing')
 require(common=={0,4,6,8,10,12,14,16,18,20,24},'common E set')
 canonical=[(E,m) for E,m in full if 0<E<12 or (E==12 and (m&127)==7)]
 canonical.sort();require(len(canonical)==12158,'canonical count')
 masks={m for E,m in canonical};require(len(masks)==12158,'duplicate canonical mask')
 require(not any((ALL^m) in masks for m in masks),'canonical complement duplicated')
 require(len(fs)==2*len(masks)+2,'full complement quotient')
 for E,m in canonical:
  for block in range(3):
   x=(m>>(7*block))&127;require(4*(x&7).bit_count()+3*(x>>3).bit_count()==E,'complex mask weight')
  x=m>>21;require(4*(x&1)+8*((x>>1)&1)+6*(x>>2).bit_count()==E,'real mask weight')
 return canonical

def disk(x):
 c=cj(x['disk']['center']);r=pf(x['disk']['radius']);require(r>0,'positive radius');return c,r
def separated(ds):
 for i,(z,r) in enumerate(ds):
  for w,t in ds[:i]:require(dis(z,r,w,t),'disk separation')
def intervals(c,r,kind):
 if kind=='linear_real':return [(-(c[0]+r),-(c[0]-r)),(F(1),F(1))]
 def square(a,b):return (F(0) if a<=0<=b else min(a*a,b*b),max(a*a,b*b))
 re=square(c[0]-r,c[0]+r);im=square(c[1]-r,c[1]+r)
 return [(re[0]+im[0],re[1]+im[1]),(-2*(c[0]+r),-2*(c[0]-r)),(F(1),F(1))]

def verify(b,L):
 require(b['format']=='H7-a3b4-mixed-rational-rouche-v1' and b['H']==H,'target H')
 scale=int(b['scale']);require(0<scale<=10**1000,'positive bounded scale')
 require(len(b['mothers'])==4 and len(b['lifts'])==25,'root record counts')
 moms={};md=[]
 for m in b['mothers']:
  k=m['block'];require(type(k)is int and 0<=k<=3 and k not in moms,'unique mother block')
  require(m['kind']==('complex_mother' if k<3 else 'real_mother'),'mother kind binding')
  z,r=disk(m)
  if k<3:require(z[1]>r,'positive nonreal mother');md.append(((z[0],-z[1]),r))
  else:require(z[1]==0 and z[0]+r<0,'negative real mother')
  lo,hi=mother(H,z,r,scale);require(lo<hi and pf(m['lhs'])==lo and pf(m['rhs'])==hi,'mother Rouche')
  moms[k]=(z,r);md.append((z,r))
 require(len(md)==7,'seven mother disks');separated(md)
 seen=set();slots=set();td=[];group_intervals={};slot_counts=Counter()
 for x in b['lifts']:
  i=x['source_group_index'];require(type(i)is int and i not in seen and 0<=i<25,'unique group index');seen.add(i)
  if i<21:k=i//7;fam='a' if i%7<3 else 'b';kind='quadratic'
  else:k=3;fam='a' if i<23 else 'b';kind='linear_real' if i==21 else 'quadratic'
  q=3 if fam=='a' else 4
  require(x['block']==k and type(x['block'])is int and x['family']==fam and x['kind']==kind and x['exponent']==q and type(x['exponent'])is int,'group mapping binding')
  j=x['lift_index'];require(type(j)is int and 0<=j<q and (k,fam,j) not in slots,'unique lift slot');slots.add((k,fam,j));slot_counts[(k,fam)]+=1
  z,r=disk(x)
  if kind=='linear_real':require(z[1]==0 and z[0]+r<0,'negative real lift')
  else:
   require(abs(z[1])>r,'nonreal quadratic lift')
   if k==3:require(z[1]>r,'real mother canonical positive lift')
   td.append(((z[0],-z[1]),r))
  lo,hi=lift(*moms[k],z,r,scale,q);require(lo<hi and pf(x['lhs'])==lo and pf(x['rhs'])==hi,'lift Rouche')
  td.append((z,r));group_intervals[i]=intervals(z,r,kind)
 require(slot_counts==Counter({(k,f):n for k in range(4) for f,n in [('a',3 if k<3 else 2),('b',4 if k<3 else 2)]}),'family root counts')
 require(len(td)==49 and len(seen)==25,'49-root coverage');separated(td)
 require(b['coverage']=={'root_disks':49,'quadratic_groups':24,'real_linear_groups':1,'pairwise_disjoint':True},'coverage metadata')
 canonical=layouts();counts=Counter(E for E,m in canonical)
 require(len(L['rows'])==5 and len({r['E'] for r in L['rows']})==5,'unique ledger rows')
 for r in L['rows']:
  E=r['E'];require(type(E)is int and E in counts and r['count']==counts[E] and type(r['count'])is int,'ledger count')
  if E<12:
   uv=[[u,v] for u in range(4) for v in range(5) if 4*u+3*v==E]
   real=[[e,n,w] for e in range(2) for n in range(2) for w in range(3) if 4*e+8*n+6*w==E]
   require(r['complex_uv']==uv and r['real_choices']==real,'ledger local choices')
  else:require(r['terminal_fixed_first_complex_all_a']is True,'terminal convention')
 require(L['total']==len(canonical) and type(L['total'])is int,'ledger total')
 require(L['groups']=={'quadratic':24,'real_linear':1},'ledger groups')
 require(L['degree_weights']=={'complex':'4u+3v','real':'4eps+8eta+6w'},'ledger weights')
 linear=group_intervals[21];require(linear[0][0]>0 and linear[1]==(1,1),'linear coefficient interval')
 enc=json.dumps(canonical,separators=(',',':')).encode()
 return {'status':'exact_frontend_and_independent_masks_passed','root_disks':49,'mother_disks':7,'layout_counts':dict(counts),'layouts':len(canonical),'canonical_E_mask_sha256':hashlib.sha256(enc).hexdigest(),'linear_constant_interval':[str(v) for v in linear[0]],'scope':'Complete necessary root-layout coverage; no coefficient exclusions.'},canonical,group_intervals

def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--ledger',type=Path,required=True);p.add_argument('--record',type=Path,required=True);p.add_argument('--masks',type=Path);a=p.parse_args()
 d,m,_=verify(json.loads(a.balls.read_text()),json.loads(a.ledger.read_text()));d.update(balls_sha256=sha(a.balls),ledger_sha256=sha(a.ledger),verifier_sha256=sha(__file__),helper_sha256=sha(Path(__file__).with_name('h7_exact_common.py')))
 if a.masks:a.masks.write_text(json.dumps(m,separators=(',',':'))+'\n')
 a.record.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
if __name__=='__main__':main()
