"""Independent full binary-pair enumeration of the d=1..9 prefix ledger."""
from pathlib import Path
from itertools import product
import argparse,json,hashlib

def verify(path):
 data=json.loads(Path(path).read_text());rows=data['rows']
 assert data['max_d']==9 and len(rows)==9
 assert all(type(r['d']) is int for r in rows)
 assert {r['d'] for r in rows}==set(range(1,10))
 counts=[]
 for row in rows:
  d=row['d'];accepted=[];survivors=[]
  for bits in product((0,1),repeat=2*(d-1)):
   a=[1]+list(bits[:d-1]);b=[1]+list(bits[d-1:])
   c=[0]*(2*d-1)
   for i,x in enumerate(a):
    for j,y in enumerate(b):c[i+j]+=x*y
   if any(x not in (0,1) for x in c[:d]):continue
   if (c[d] if d<len(c) else 0)!=0:continue
   assert all(a[i]+b[i]<=1 for i in range(1,d))
   code=''.join(str(a[i]+2*b[i]) for i in range(1,d));accepted.append(code)
   bad=max(c)>1
   for i in range(1,d):
    baseline=c[d+i] if d+i<len(c) else 0
    if baseline>=1 and a[i]+b[i]>0:bad=True
   if not bad:survivors.append(code)
  accepted.sort();survivors.sort();assert len(set(accepted))==len(accepted)
  assert accepted==row['valid_codes'] and len(accepted)==row['valid_count']
  assert row['raw_prefix_codes']==3**(d-1)
  assert hashlib.sha256(json.dumps(accepted,separators=(',',':')).encode()).hexdigest()==row['code_sha256']
  f=row['strict_necessary_filter'];assert survivors==f['survivor_codes'] and len(survivors)==f['survivor_count']
  assert hashlib.sha256(json.dumps(survivors,separators=(',',':')).encode()).hexdigest()==f['survivor_code_sha256']
  counts.append([d,4**(d-1),len(accepted),len(survivors)])
 return {'status':'all_binary_pair_prefixes_verified','range':[1,9],'counts_d_raw4_legal_survivor':counts,'ledger_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Necessary low-end patterns and two necessary global budget filters; not existence of extensions.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args();v=verify(a.input);a.record.write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v))
