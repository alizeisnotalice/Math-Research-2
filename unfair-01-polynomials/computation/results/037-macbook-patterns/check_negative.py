"""Exercise the production stream verifier on a complete tiny valid ledger."""
from pathlib import Path
import copy,gzip,json,subprocess,sys,tempfile,os,hashlib
assert os.path.ismount('/Volumes/data')
root=Path(__file__).resolve().parent
temp=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round037-negative-portable');temp.mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(dir=temp) as work:
    work=Path(work);p=work/'valid.json.gz'
    subprocess.run([sys.executable,str(root/'generate_one.py'),'--p','2','--q','5','--out',str(p)],check=True,capture_output=True)
    d=json.loads(gzip.decompress(p.read_bytes()))
    def run(x,corrupt=False):
        b=gzip.compress(json.dumps(x,separators=(',',':')).encode(),mtime=0)
        if corrupt:b=b[:-4]+b'BAD!'
        return subprocess.run([sys.executable,str(root/'verify_pattern.py'),'2','5'],input=b,capture_output=True)
    baseline=run(d);assert baseline.returncode==0,baseline.stderr
    mutators={'drop_row':lambda x:x['rows'].pop(),'duplicate_row':lambda x:x['rows'].append(copy.deepcopy(x['rows'][0])),
      'false_units':lambda x:x['rows'][0]['unit_indices'].pop(),'false_gap':lambda x:x['rows'][0]['trace'].insert(0,['zero-gap',0,1,0,1,100]),
      'wrong_count':lambda x:x.__setitem__('unknown',x['unknown']+1),'wrong_parameter':lambda x:x['rows'][0].__setitem__('p',3)}
    rejected={}
    for name,mut in mutators.items():
        x=copy.deepcopy(d);mut(x);r=run(x);assert r.returncode;rejected[name]=r.stderr.decode().splitlines()[-1]
    r=run(d,True);assert r.returncode;rejected['corrupt_gzip']=r.stderr.decode().splitlines()[-1]
    print(json.dumps({'status':'passed','positive':json.loads(baseline.stdout.decode().splitlines()[-1]),'rejected':rejected,'verifier_sha256':hashlib.sha256((root/'verify_pattern.py').read_bytes()).hexdigest(),'python':sys.version},indent=2))
