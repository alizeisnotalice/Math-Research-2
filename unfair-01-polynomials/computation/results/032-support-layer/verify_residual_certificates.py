"""Recheck saved strict diagonal injections; do not run the certificate searcher."""
import hashlib, importlib.util, json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main():
    path=HERE.parents[1]/'diagonal_domination.py'
    spec=importlib.util.spec_from_file_location('strict_diagonal_verifier',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    data=json.loads((HERE/'residual-46-mapping.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==data['source_sha256']['diagonal_domination.py']
    seen=set();proved=[];unknown=[]
    for row in data['minimal_unknown_inputs']:
        i,j=row['I'],row['J'];key=(tuple(i),tuple(j))
        assert key not in seen and i==sorted(set(i)) and j==sorted(set(j)) and i[0]==j[0]==0
        seen.add(key);cert=row['direct_diagonal_domination']
        if cert is None:unknown.append({'I':i,'J':j})
        else:
            assert mod.verify_certificate(i,j,cert)
            proved.append({'I':i,'J':j})
    assert len(seen)==46 and len(proved)==39 and len(unknown)==7
    out={'status':'39_saved_strict_certificates_independently_recomputed','proved_count':len(proved),
         'unresolved_count':len(unknown),'unresolved':unknown,'searcher_called':False,
         'verifier_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'scope':'Only the 39 saved certificates are accepted; the 7 absent certificates remain unknown.'}
    (HERE/'parent-residual-verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print({'proved':len(proved),'unknown':len(unknown)})
if __name__=='__main__':main()
