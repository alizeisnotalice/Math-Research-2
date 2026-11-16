"""Bind every one of the 46 residual supports to a rechecked exact certificate."""
import argparse,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(x):return tuple(x['I']),tuple(x['J'])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--last-certificate',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args()
    old=a.topic/'computation/results/032-support-layer';mapping=old/'residual-46-mapping.json';data=json.loads(mapping.read_text())
    diagonal=a.topic/'computation/diagonal_domination.py'
    assert sha(diagonal)==data['source_sha256']['diagonal_domination.py']
    spec=importlib.util.spec_from_file_location('saved_diagonal_verifier',diagonal);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    remaining={};rows=[];seen=set()
    for row in data['minimal_unknown_inputs']:
        k=key(row);assert k not in seen;seen.add(k)
        cert=row['direct_diagonal_domination']
        if cert is None:remaining[k]=row
        else:
            assert mod.verify_certificate(row['I'],row['J'],cert)
            rows.append({'I':row['I'],'J':row['J'],'method':'strict_diagonal_injection','source_unknown':row['source_unknown']})
    assert len(seen)==46 and len(rows)==39 and len(remaining)==7
    verifier=a.topic/'computation/verify_positive_certificate.py'
    for path in sorted((old/'positive').glob('*.json'))+[a.last_certificate]:
        d=json.loads(path.read_text());k=key(d['case']);assert k in remaining
        row=remaining.pop(k)
        run=subprocess.run([sys.executable,str(verifier),str(path)],capture_output=True,text=True,check=True)
        result=json.loads(run.stdout);assert result['status']=='verified'
        rows.append({'I':row['I'],'J':row['J'],'method':'positive_polynomial_identity','source_unknown':row['source_unknown'],'certificate_path':str(path),'certificate_sha256':sha(path),'verification':result})
    assert not remaining and len(rows)==46
    result={'status':'all_46_residual_supports_exactly_certified','direct_count':39,'positive_count':7,'unknown':0,'mapping_sha256':sha(mapping),'verifier_sha256':sha(verifier),'diagonal_sha256':sha(diagonal),'rows':rows,'scope':'Residual certificates only; finite layer conclusion also requires full support coverage and source-bound verification receipts.'}
    a.record.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','direct_count','positive_count','unknown']}))
if __name__=='__main__':main()
