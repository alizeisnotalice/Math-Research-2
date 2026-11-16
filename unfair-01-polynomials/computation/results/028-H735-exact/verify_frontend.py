"""Independent exact H7(3,5) root cover and 29-bit necessary layout cover."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse,json,hashlib,gzip
from exact_complex import cj,pf,mother,lift,dis
H=[1,1,1,0,0,0,1,1];ALL=(1<<29)-1

def req(x,s):
 if not x:raise ValueError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def layouts():
 local=[]
 for block in range(3):
  d={}
  for bits in range(256):d.setdefault(5*(bits&7).bit_count()+3*(bits>>3).bit_count(),[]).append(bits<<(8*block))
  local.append(d)
 d={}
 for bits in range(32):d.setdefault(5*(bits&1)+10*((bits>>1)&1)+3*((bits>>2)&1)+6*(bits>>3).bit_count(),[]).append(bits<<24)
 local.append(d);common=set.intersection(*(set(x) for x in local));full=[]
 for E in sorted(common):full.extend((E,a|b|c|r) for a,b,c,r in product(*(x[E] for x in local)))
 fs=set(full);req(len(fs)==len(full)==238126,'full unique masks');req(all((30-E,ALL^m) in fs for E,m in fs),'complement closure')
 canon=sorted((E,m) for E,m in full if 0<E<15 or (E==15 and m&255==7))
 req(len(canon)==119062 and len({m for E,m in canon})==119062,'canonical mask count')
 cs=set(canon);req(all((30-E,ALL^m) not in cs for E,m in cs),'canonical quotient');req(len(fs)==2*len(cs)+2,'empty/full separated')
 for E,m in canon:
  for k in range(3):
   x=(m>>(8*k))&255;req(5*(x&7).bit_count()+3*(x>>3).bit_count()==E,'complex weight')
  x=m>>24;req(5*(x&1)+10*((x>>1)&1)+3*((x>>2)&1)+6*(x>>3).bit_count()==E,'real weight')
 req(dict(Counter(E for E,m in canon))=={3:125,5:27,6:2000,8:3375,9:2000,10:27,11:54000,12:125,13:3375,14:54000,15:8},'E counts')
 return canon

def disk(x):
 z=cj(x['disk']['center']);r=pf(x['disk']['radius']);req(r>0,'positive disk radius');return z,r
def sep(ds):
 for i,(z,r) in enumerate(ds):
  for c,t in ds[:i]:req(dis(z,r,c,t),'root disk separation')
def box(z,r,kind):
 if kind=='linear':return [(-z[0]-r,-z[0]+r),(F(1),F(1))]
 def sq(l,h):return (F(0) if l<=0<=h else min(l*l,h*h),max(l*l,h*h))
 re=sq(z[0]-r,z[0]+r);im=sq(z[1]-r,z[1]+r)
 return [(re[0]+im[0],re[1]+im[1]),(-2*(z[0]+r),-2*(z[0]-r)),(F(1),F(1))]

def verify(b):
 req(b['format']=='H7-a3b5-mixed-rational-rouche-v1' and b['H']==H and b['a']==3 and b['b']==5,'fixed target');scale=int(b['scale']);req(0<scale<=10**1000,'scale')
 support=[i for i,x in enumerate(H) if x];C=[0]*57
 for i in support:
  for j in support:C[3*i+5*j]+=1
 req(C[0]==C[-1]==1 and set(C)=={0,1} and sum(C)==25,'binary target')
 req(len(b['mothers'])==4 and len(b['groups'])==29,'source counts');moms={};md=[]
 for m in b['mothers']:
  k=m['block'];req(type(k)is int and 0<=k<=3 and k not in moms,'unique mother block');z,r=disk(m)
  if k<3:req(z[1]>r,'positive nonreal mother');md.append(((z[0],-z[1]),r))
  else:req(z[1]==0 and z[0]+r<0,'negative real mother')
  L,R=mother(H,z,r,scale);req(L<R and pf(m['lhs'])==L and pf(m['rhs'])==R,'mother Rouche');md.append((z,r));moms[k]=(z,r)
 req(len(md)==7,'seven mother roots');sep(md);seen=set();td=[];f={};counts=Counter();families=Counter();standard=0
 for g in b['groups']:
  i=g['index'];req(type(i)is int and i not in seen and 0<=i<29,'unique source group');seen.add(i)
  if i<24:k=i//8;off=i%8;fam='a' if off<3 else 'b';slot=off if fam=='a' else off-3;kind='quadratic'
  else:
   k=3;fam='a' if i<26 else 'b';slot=i-24 if fam=='a' else i-26;kind='linear' if i in (24,26) else 'quadratic'
  req(type(g['block'])is int and g['block']==k and g['family']==fam and type(g['slot'])is int and g['slot']==slot and g['kind']==kind,'group binding')
  q=3 if fam=='a' else 5;z,r=disk(g)
  if kind=='linear':req(z[1]==0 and z[0]+r<0,'negative real lift')
  else:
   req(abs(z[1])>r,'nonreal lift')
   if k==3:req(z[1]>r,'real mother positive representative')
   td.append(((z[0],-z[1]),r))
  L,R=lift(*moms[k],z,r,scale,q);req(L<R and pf(g['lhs'])==L and pf(g['rhs'])==R,'lift Rouche');td.append((z,r));f[i]=box(z,r,kind);counts[k,fam]+=1;families[fam]+=len(f[i])-1
  if fam=='a':standard|=1<<i
 req(counts==Counter({(k,fam):(3 if fam=='a' else 5) if k<3 else (2 if fam=='a' else 3) for k in range(4) for fam in ('a','b')}),'block counts')
 req(families=={'a':21,'b':35} and len(td)==56,'root family degrees');sep(td)
 req(b['coverage']=={'root_disks':56,'quadratic_groups':27,'real_linear_groups':2,'pairwise_disjoint':True},'coverage metadata')
 for i in (24,26):req(f[i][0][0]>0 and f[i][-1]==(1,1),'linear sign')
 rows=layouts();req((15,standard) in rows,'standard included');encoded=json.dumps(rows,separators=(',',':')).encode()
 return {'status':'exact_root_and_mask_frontend_passed','degree':56,'binary_terms':25,'mother_roots':7,'target_roots':56,'quadratic_groups':27,'linear_groups':2,'layouts':len(rows),'counts_by_E':dict(Counter(E for E,m in rows)),'canonical_sha256':hashlib.sha256(encoded).hexdigest(),'standard_mask':standard,'scope':'Necessary coverage only; no coefficient exclusions.'},rows,f

def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--record',type=Path,required=True);p.add_argument('--masks',type=Path);p.add_argument('--source-ledger',type=Path);a=p.parse_args();v,rows,_=verify(json.loads(a.balls.read_text()))
 if a.source_ledger:
  opener=gzip.open if a.source_ledger.suffix=='.gz' else open
  with opener(a.source_ledger,'rt') as f:source=[json.loads(line) for line in f]
  # Source schema is inspected separately; only the independently regenerated
  # complete set may establish coverage.
  pairs=[]
  for x in source:
   bits=x['mask'];req(type(bits)is str and len(bits)==29 and set(bits)<=set('01'),'external bit string')
   ids=[i for i,ch in enumerate(bits) if ch=='1'];req(x['selected_group_indices']==ids,'external selected ids');req(type(x['E'])is int,'external integer E');pairs.append((x['E'],sum(1<<i for i in ids)))
  req(len(pairs)==len(set(pairs))==len(rows) and set(pairs)==set(rows),'external ledger full set');v['source_ledger_sha256']=sha(a.source_ledger)
 v.update(balls_sha256=sha(a.balls),verifier_sha256=sha(__file__),helper_sha256=sha(Path(__file__).with_name('exact_complex.py')))
 if a.masks:a.masks.write_text(json.dumps(rows,separators=(',',':'))+'\n')
 a.record.write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v))
if __name__=='__main__':main()
