"""Certify all angular cells outside the one with the unique triple root."""
import json,hashlib
from pathlib import Path
from flint import arb,ctx
ctx.prec=160;p=Path(__file__).resolve().parent;load=lambda n:json.loads((p/n).read_text());c,f,m=map(load,['certificate.json','fine-certificate.json','micro-certificate.json']);cs,fs,ms=map(load,['scan.json','fine-scan.json','micro-scan.json']);sel=load('selected.json');mini=load('micro-selected.json');d=load('derivative-cells.json');cd=load('coarse-derivative-cells.json');md=load('micro-derivative-cells.json');tr=load('triple.json');loc=load('local-exact.json');upper=arb(tr['area_upper']).upper()
for rec,name in [(c,'scan.json'),(f,'fine-scan.json'),(m,'micro-scan.json')]:assert rec['source_scan_sha256']==hashlib.sha256((p/name).read_bytes()).hexdigest() and not rec['failed']
assert [x['i'] for x in cs['results']]==list(range(10000)) and [x['i'] for x in fs['results']]==[100*i+j for i in sel for j in range(100)] and [x['i'] for x in ms['results']]==[10*i+j for i in mini for j in range(10)]
assert len(set(sel))==len(sel) and len(set(mini))==len(mini) and set(mini)<=set(x['i'] for x in fs['results'])
D={x['i']:arb(x['lower']).lower() for x in d['results']};C={x['i']:arb(x['lower']).lower() for x in cd['results']};T={x['i']:arb(x['lower']).lower() for x in md['results']};assert len(D)==len(d['results']) and len(C)==len(cd['results']) and len(T)==len(md['results']) and set(T)==set(x['i'] for x in ms['results'])
assert loc['interval_cell']==344314 and loc['contacts']==[0,7,26] and arb(loc['proof'][0]['derivative'])<0 and arb(loc['proof'][1]['derivative'])>0 and arb(loc['hinge_nonzero_margin_lower'])>0
bounds=[]
for i,x in enumerate(c['cell_lowers']):
 if i not in sel:bounds.append((C[i] if i in C else (arb(x).lower()-arb(c['angle_error'])).lower(),'coarse',i))
for z,x in zip(fs['results'],f['cell_lowers']):
 i=z['i']
 if i not in mini:bounds.append((D[i] if i in D else (arb(x).lower()-arb(f['angle_error'])).lower(),'fine',i))
for i,x in T.items():
 if i!=loc['interval_cell']:bounds.append((x,'micro',i))
assert 344314 in T and 34431 in mini
lower,layer,index=min(bounds,key=lambda x:x[0]);assert lower>upper,(layer,index,lower,upper)
out={'exact_global_minimum':'triple.json: quadratic and isolated tan-half-angle root; area at three-contact pose','global_exclusion_margin_lower':str((lower-upper).lower()),'outer_bottleneck':{'layer':layer,'index':index},'three_contact_angles_covered_by_sign_definite_balanced_duals':True,'both_hands_and_full_angle':True,'precision_bits':ctx.prec};(p/'joined-exact.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
