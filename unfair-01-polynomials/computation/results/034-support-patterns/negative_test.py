"""Reject a support omission even after the advertised counts are changed."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert Path('/Volumes/data').is_mount();a.out.mkdir(parents=True,exist_ok=True)
    d=json.loads(a.input.read_text());r=next(x for x in d['rows'] if (x['p'],x['q'])==(2,5))
    r['no_mod_p_singleton_survivors'].pop();r['no_mod_p_singleton_count']-=1;d['summary']['no_singleton_survivors']-=1
    bad=a.out/'missing-support.json';bad.write_text(json.dumps(d)+'\n')
    verifier=Path(__file__).with_name('audit_short_g_singletons.py')
    run=subprocess.run([sys.executable,str(verifier),'--input',str(bad),'--record',str(a.out/'invalid-receipt.json')],text=True,capture_output=True)
    assert run.returncode!=0 and 'concrete survivor omission/addition' in run.stderr
    out={'status':'rejected','mutation':'Delete one p2q5 survivor with synchronized row/top count decrements','returncode':run.returncode,'stderr_tail':run.stderr[-500:],'verifier_sha256':hashlib.sha256(verifier.read_bytes()).hexdigest()}
    (a.out/'negative-test-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'])
if __name__=='__main__':main()
