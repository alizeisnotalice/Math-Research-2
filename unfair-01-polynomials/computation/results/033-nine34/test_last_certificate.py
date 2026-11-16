"""Production verifier regressions against the new exact positive identity."""
import argparse,copy,hashlib,json,subprocess,sys
from fractions import Fraction
from pathlib import Path
def main():
    p=argparse.ArgumentParser();p.add_argument('--verifier',type=Path,required=True);p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert Path('/Volumes/data').is_mount();a.out.mkdir(parents=True,exist_ok=True)
    base=json.loads(a.certificate.read_text());tests={}
    def execute(name,data,accept):
        path=a.out/(name+'.json');path.write_text(json.dumps(data)+'\n')
        r=subprocess.run([sys.executable,str(a.verifier),str(path)],capture_output=True,text=True)
        assert (r.returncode==0)==accept,(name,r.stdout,r.stderr)
        tests[name]={'returncode':r.returncode,'expected_accept':accept,'stderr':r.stderr.strip()}
    execute('valid',base,True)
    d=copy.deepcopy(base);term=d['certificate']['multipliers'][0];q=Fraction(*term['coefficient'])+Fraction(1,10**30);term['coefficient']=[q.numerator,q.denominator]
    execute('tiny_exact_identity_perturbation',d,False)
    d=copy.deepcopy(base);d['certificate']['positive_atoms'][0]['coefficient'][0]*=-1
    execute('negative_positive_atom',d,False)
    d=copy.deepcopy(base);d['strict_variables'].pop()
    execute('missing_strict_variable',d,False)
    out={'status':'passed','certificate_sha256':hashlib.sha256(a.certificate.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(a.verifier.read_bytes()).hexdigest(),'tests':tests}
    (a.out/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
