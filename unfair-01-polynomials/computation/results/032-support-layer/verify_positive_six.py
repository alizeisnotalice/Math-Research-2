"""Bind six positive identities to distinct residual supports and run exact verifier."""
import hashlib,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def key(row):return tuple(row['I']),tuple(row['J'])
def main():
    remaining=json.loads((HERE/'parent-residual-verification.json').read_text())['unresolved']
    expected={key(r) for r in remaining};seen=set();receipts=[]
    verifier=HERE.parents[1]/'verify_positive_certificate.py'
    for path in sorted((HERE/'positive').glob('*.json')):
        data=json.loads(path.read_text());k=key(data['case'])
        assert k in expected and k not in seen;seen.add(k)
        run=subprocess.run([sys.executable,str(verifier),str(path)],capture_output=True,text=True,check=True)
        receipt=json.loads(run.stdout);assert receipt['status']=='verified'
        receipts.append({'path':str(path.relative_to(HERE)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                         'verification':receipt})
    assert len(seen)==6 and len(expected-seen)==1
    out={'status':'six_new_positive_identities_exact_verified','direct_already_verified':39,
         'positive_verified':6,'total_accepted':45,'remaining':[r for r in remaining if key(r) not in seen],
         'verifier_sha256':hashlib.sha256(verifier.read_bytes()).hexdigest(),'receipts':receipts,
         'cap9_total34_closed':False}
    (HERE/'parent-positive-six.json').write_text(json.dumps(out,indent=2)+'\n')
    print({'accepted':45,'remaining':out['remaining']})
if __name__=='__main__':main()
