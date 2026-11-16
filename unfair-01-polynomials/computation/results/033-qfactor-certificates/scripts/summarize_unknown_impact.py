#!/usr/bin/env python3
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
def k(a):return hashlib.sha256(json.dumps(a,separators=(',',':')).encode()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--search-out',type=Path,required=True);p.add_argument('--certificates',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args();c=json.loads(a.certificates.read_text());unknown={r['factor_sha256'] for r in c['certificates'] if r['certificate']['kind'].startswith('unknown')};tot=Counter();ids=[]
 with (a.search_out/'h-coverage.jsonl').open() as f:
  for line in f:
   r=json.loads(line)
   if any(k(q) in unknown for q in r['factor_coefficients_low_to_high']):
    ids.append(r['h_id']);tot['H_touching_unknown']+=1
    for name,value in r['allocation_counts'].items():tot['allocations_'+name]+=value
 out={'status':'unknown_impact_counted','unknown_factor_count':len(unknown),'H_touching_unknown_count':len(ids),'H_ids_sha256':hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest(),'H_ids':ids,'allocation_counts_for_H_touching_unknown':dict(tot),'claim_boundary':'These H/allocation rows lack complete irreducibility certification because at least one listed factor remains unknown; they are not excluded by the certified subset.'};a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({k:out[k] for k in ('unknown_factor_count','H_touching_unknown_count','H_ids_sha256')},sort_keys=True))
if __name__=='__main__':main()
