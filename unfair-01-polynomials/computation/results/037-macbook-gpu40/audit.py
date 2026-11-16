from pathlib import Path
import json,hashlib,time,numpy as np
r=Path(__file__).resolve().parent;data=json.loads((r/'inputs.json').read_text());events=[]
for l in (r/'stdout.jsonl').read_text().splitlines():
 if l.startswith('{'):events.append(json.loads(l))
assert events[-1]['event']=='complete'
root_product_error=0.
for t in data['tasks']:
 assert len(t['C'])==41 and t['C'][0]==t['C'][-1]==1 and set(t['C'])<={0,1}
 assert t['groups']==len(t['roots']) and t['masks']==2**(t['groups']-1)
 allroots=[complex(*z) for group in t['roots'] for z in group];assert len(allroots)==40
 full=np.polynomial.polynomial.polyfromroots(allroots);full/=full[0]
 err=float(np.max(abs(full-np.array(t['C']))));root_product_error=max(root_product_error,err);assert err<1e-6
cursors={t['id']:0 for t in data['tasks']};count=0;maxerr=0;best=1e100;nsamples=0;out=[];product_error=0.
for e in events:
 if e['event']!='batch':continue
 t=data['tasks'][e['task']];assert e['first']==cursors[t['id']] and e['count']==e['last']-e['first']+1 and e['gpu_arrays_verified'];cursors[t['id']]=e['last']+1;count+=e['count']
 for c in e['cases']:
  bits=[0]+list(map(int,format(c['mask'],f"0{t['groups']-1}b"))); aa=np.array([1.]);bb=np.array([1.])
  for bit,group in zip(bits,t['roots']):
   roots=[complex(*v) for v in group];poly=np.polynomial.polynomial.polyfromroots(roots);poly/=poly[0];assert np.max(abs(poly.imag))<1e-8;poly=poly.real
   if bit:aa=np.polynomial.polynomial.polymul(aa,poly)
   else:bb=np.polynomial.polynomial.polymul(bb,poly)
  aa=np.pad(aa,(0,41-len(aa)));bb=np.pad(bb,(0,41-len(bb)));err=max(np.max(abs(aa-c['A'])),np.max(abs(bb-c['B'])));maxerr=max(maxerr,float(err));assert err<1e-4;product=np.polynomial.polynomial.polymul(aa,bb);target=np.pad(t['C'],(0,len(product)-len(t['C'])));product_error=max(product_error,float(np.max(abs(product-target))));assert np.max(abs(product-target))<1e-5;score=max(0.,-min(aa),-min(bb),max(aa)-1,max(bb)-1);best=min(best,float(score));nsamples+=1
  if score<1e-4:out.append({'task':t['id'],'mask':c['mask'],'score':float(score),'gpu_score':c['score'],'coefficient_error':float(err),'A':aa.tolist(),'B':bb.tolist()})
assert all(cursors[t['id']]==t['masks'] for t in data['tasks'])
report={'status':'numerical_samples_reconstructed','masks':count,'samples':nsamples,'max_gpu_cpu_coefficient_error':maxerr,'max_root_product_error':root_product_error,'max_sample_product_error':product_error,'minimum_cpu_score':best,'near_candidates':out,'scope':'float64 reconstruction from saved approximate roots only; no exact exclusion','numpy':np.__version__,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(r/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='near_candidates'}));print('near',len(out))
