#!/usr/bin/env python3
import argparse,copy,json,subprocess,sys
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument('--verifier',type=Path,required=True);p.add_argument('--certificates',type=Path,required=True);p.add_argument('--search-out',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();base=json.loads(a.certificates.read_text());rabin=next(r for r in base['certificates'] if r['certificate']['kind']=='irreducible_mod_prime_rabin');tests={}
 for name,mutate in {
  'deleted_gcd':lambda d,r:r['certificate']['rabin_gcd_checks'].pop(),
  'composite_p':lambda d,r:r['certificate'].__setitem__('prime',4),
  'wrong_degree':lambda d,r:r['certificate'].__setitem__('degree',r['certificate']['degree']+1),
  'duplicate_certificate_row':lambda d,r:d['certificates'].append(copy.deepcopy(d['certificates'][0]))
 }.items():
  d=copy.deepcopy(base);rr=next(x for x in d['certificates'] if x['factor_sha256']==rabin['factor_sha256']);mutate(d,rr);path=a.out/(name+'.json');path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(d,separators=(',',':'))+'\n');record=a.out/(name+'-receipt.json');q=subprocess.run([sys.executable,str(a.verifier),'--search-out',str(a.search_out),'--certificates',str(path),'--record',str(record)],text=True,capture_output=True)
  tests[name]={'returncode':q.returncode,'rejected':q.returncode!=0,'stderr_tail':q.stderr[-400:],'stdout_tail':q.stdout[-400:]}
 if not all(x['rejected'] for x in tests.values()):raise SystemExit('a mutation was accepted')
 (a.out/'negative-tests.json').write_text(json.dumps({'status':'all_negative_tests_rejected','tests':tests},indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'all_negative_tests_rejected','tests':tests},sort_keys=True))
if __name__=='__main__':main()
