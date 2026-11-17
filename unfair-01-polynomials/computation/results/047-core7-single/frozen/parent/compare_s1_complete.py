#!/usr/bin/env python3
import argparse,hashlib,json,math,gzip
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--records',type=Path,required=True);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 expected=sorted((T,(T*T-D)//4,D) for T in range(-16,17) for D in range(1,(17-abs(T))**2) if (T*T-D)%4==0 and math.isqrt(D)**2!=D)
 expected.sort(key=lambda t:(t[0],t[2]));assert len(expected)==680
 parent={}
 for line in a.parent.read_text().splitlines():
  r=json.loads(line);k=tuple(r['TN']);assert k not in parent;parent[k]=r
 assert set(parent)==set(expected);source_hash=sha(a.source);bindings=[];counts={k:0 for k in ['c2','c3','c4','prefix4','c5','prefix6']};sole=[]
 assert len(list(a.records.glob('*.json')))==680
 for ix,triple in enumerate(expected):
  p=a.records/f'{ix:04}.json';x=json.loads(p.read_text());assert x['index']==ix and tuple(x['TNDelta'])==triple and x['returncode']==0 and x['status']=='completed' and x['source_sha256']==source_hash
  y=x['result'];assert y['tn_start']==ix and y['tn_count']==1 and y['input_TN_total']==680 and y['script_sha256']==source_hash and y['counts']['TN']==1
  z=parent[triple]
  for k in counts:assert y['counts'][k]==z['counts'][k];counts[k]+=y['counts'][k]
  states=sorted([v['c2'],v['c3'],v['c4'],v['c5'],v['G'],v['Q']] for v in y['survivors']);assert states==sorted(z['passed6'])
  if states:sole.append({'TN':triple,'states':states})
  bindings.append({'index':ix,'TN':triple,'sha256':sha(p)})
 out={'status':'passed_all_680_input_triples_stage_counts_and_survivors','counts':counts,'sole':sole,'parent_ledger_sha256':sha(a.parent),'search_source_sha256':source_hash,'record_bindings':bindings};a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','counts','sole']}))
if __name__=='__main__':main()
