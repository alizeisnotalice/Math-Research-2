#!/usr/bin/env python3
"""Exact derivative-period audit of the 20 round005/006 dim-15 certificates."""
import argparse,hashlib,json,sys,time
P=64;FULL=(1<<64)-1
def rank(vs,n=64):
 a=vs[:];r=0
 for b in range(n-1,-1,-1):
  p=next((i for i in range(r,len(a)) if a[i]>>b&1),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r]
  for i in range(len(a)):
   if i!=r and a[i]>>b&1:a[i]^=a[r]
  r+=1
 return r
def basis(vs):
 out=[]
 for x in vs:
  if rank(out+[x],6)>len(out):out.append(x)
 return out
def fwht(a):
 a=a[:];h=1
 while h<len(a):
  for p in range(0,len(a),2*h):
   for i in range(p,p+h):a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
  h*=2
 return a
def seed(row):return [(bin((x&7)&row['permutation'][x>>3]).count('1')^((row['h_bits']>>(x>>3))&1))&1 for x in range(P)]
def cspace():return [FULL]+[sum(((x>>j)&1)<<x for x in range(P)) for j in range(6)]
def diff_const(mask,a):
 d=((mask>>0)&1)^((mask>>a)&1)
 return all((((mask>>x)&1)^((mask>>(x^a))&1))==d for x in range(P))
def certify(row):
 G=cspace()+row['relative_generators'];assert rank(G,64)==15
 periods=[a for a in range(P) if all(diff_const(g,a) for g in G)]
 # Independent graph-period form: labels are all generator evaluations;
 # translation by a is a graph period iff label(x+a)+label(x) is constant.
 labels=[sum(((g>>x)&1)<<j for j,g in enumerate(G)) for x in range(P)]
 graphperiod=[a for a in range(P) if len({labels[x]^labels[x^a] for x in range(P)})==1]
 assert periods==graphperiod
 out={'source_file':row['source_file'],'source_sha256':row['source_sha256'],'row_index':row['row_index'],'seed_index':row['seed_index'],'seed_permutation':row['permutation'],'seed_h_bits':row['h_bits'],'V_basis_masks':G,'P_elements':periods,'P_dimension':rank(periods,6),'P_basis':basis(periods),'graph_period_elements':graphperiod}
 if periods==[0]:
  q=seed(row);ok=True
  for code in range(1<<15):
   vals=[-1 if q[x]^sum(((code>>j)&1)*((G[j]>>x)&1) for j in range(15))&1 else 1 for x in range(P)]
   if any(abs(z)!=8 for z in fwht(vals)):ok=False;break
  assert ok
  out.update({'all_32768_slices_independently_bent':True,'graph_support_64_points_no_nonzero_translation_period':True})
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--audit',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();t=time.perf_counter();raw=open(a.audit,'rb').read();d=json.loads(raw);rows=[certify(r) for r in d['rows']];assert all(r.get('all_32768_slices_independently_bent') is True for r in rows if r['P_dimension']==0);o={'status':'passed','audit_input_sha256':hashlib.sha256(raw).hexdigest(),'certificate_count':len(rows),'rows':rows,'source_sha256':hashlib.sha256(open(__file__,'rb').read()).hexdigest(),'python':sys.version,'elapsed_seconds':time.perf_counter()-t};open(a.out,'w').write(json.dumps(o,indent=2)+'\n');print(json.dumps({'status':'passed','P_dimensions':[r['P_dimension'] for r in rows],'P0_count':sum(r['P_dimension']==0 for r in rows),'elapsed_seconds':o['elapsed_seconds']}))
if __name__=='__main__':main()
