"""Extract only the five H10 upper-half mother disks from the old 2/3 input."""
import argparse, hashlib, json
from pathlib import Path

H10=[1,1,1,0,0,0,0,0,0,1,1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 raw=a.source.read_bytes();d=json.loads(raw)
 assert d['format']=='H10-a2b3-rational-rouche-balls-v1' and d['H']==H10 and len(d['mothers'])==5
 out={'format':'H10-mother-balls-from-026-v1','H':H10,'scale':d['scale'],
      'source_root_balls_sha256':hashlib.sha256(raw).hexdigest(),
      'mothers':sorted(d['mothers'],key=lambda x:x['block']),
      'coverage':{'upper_half_mother_disks':5,'conjugate_root_disks':10,
                  'old_lifts_intentionally_omitted':True}}
 a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print(hashlib.sha256(a.output.read_bytes()).hexdigest())
if __name__=='__main__':main()
