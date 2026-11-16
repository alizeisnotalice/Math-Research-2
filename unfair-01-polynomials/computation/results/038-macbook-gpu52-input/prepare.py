from pathlib import Path
import json,os,numpy as np
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent; rng=np.random.default_rng(2026111616); tasks=[];L=64
for ti in range(24):
 c=rng.integers(0,2,53); c[0]=c[-1]=1
 if ti==0:
  c[:]=0
  for i in [0,1,3]:
   for j in [0,5,11,20,49]: c[i+j]=1
 roots=np.polynomial.polynomial.polyroots(c); unused=set(range(52)); groups=[]
 while unused:
  i=min(unused);unused.remove(i); z=roots[i]
  if abs(z.imag)<1e-8: group=[complex(z.real,0)]
  else:
   j=min(unused,key=lambda j:abs(roots[j]-z.conjugate())); assert abs(roots[j]-z.conjugate())<1e-8;unused.remove(j); group=[z,roots[j]]
  groups.append(group)
 points=1.031*np.exp(2j*np.pi*np.arange(L)/L); logs=[]
 for group in groups:
  vals=np.ones(L,dtype=complex)
  for z in group: vals*= (points-z)/(-z)
  assert np.min(abs(vals))>1e-10
  logs.append(np.log(vals))
 mat=np.array(logs)
 tasks.append({'id':ti,'C':c.tolist(),'roots':[[[z.real,z.imag] for z in g] for g in groups],'LR':mat.real.tolist(),'LI':mat.imag.tolist(),'groups':len(groups),'masks':2**(len(groups)-1)})
phase=2*np.pi*np.outer(np.arange(L),np.arange(53))/L
data={'seed':2026111616,'L':L,'batch':32768,'tasks':tasks,'DC':(np.cos(phase)/L/(1.031**np.arange(53))).tolist(),'DS':(np.sin(phase)/L/(1.031**np.arange(53))).tolist(),'scope':'numerical real root-group assignments; no exclusion or exact certificate'}
(r/'inputs.json').write_text(json.dumps(data,separators=(',',':')))
print([(t['id'],t['groups'],t['masks']) for t in tasks])
