"""Recompute every lattice pair and all fourth-row contradictions."""
import argparse,gzip,hashlib,importlib.util,json,tempfile,time,sys
from collections import Counter
from pathlib import Path
def mod(p,n):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();r=a.package;begin=time.monotonic()
    manifest=json.loads((r/'manifest.json').read_text())
    for z in manifest['files']:
        b=(r/z['path']).read_bytes();assert len(b)==z['bytes'] and sha(b)==z['sha256']
    # Temporary extraction stays next to the caller's explicit output path.
    with tempfile.TemporaryDirectory(prefix='s12-replay-',dir=a.out.parent) as tmp:
        p=Path(tmp)/'low.json';p.write_bytes(gzip.decompress((r/'data/low-ledger.json.gz').read_bytes()))
        low=mod(r/'frozen/independent_low.py','low');cd=mod(r/'frozen/independent_cd.py','cd')
        report,rows=low.main(p);ledger=json.loads(p.read_text());assert report['counts']['raw']==396738 and len(rows)==len(ledger['rows'])==7245
        counts=Counter();unit=[];survivors=[]
        for row in ledger['rows']:
            c,s,u=cd.one(row);counts.update(c);unit.extend(u);survivors.extend(s)
        assert not survivors and counts['c_raw']==760725 and counts['c_congruence']==716608 and counts['s3_pass']==19728
        assert len(unit)==counts['d_interval_candidates']==194 and all(z['C4']==2 for z in unit)
        old=json.loads(gzip.decompress((r/'data/parent-cd.json.gz').read_bytes()))
        key=lambda z:(z['prefix']['D'],z['prefix']['a_m'],z['prefix']['a_n'],z['prefix']['b_m'],z['prefix']['b_n'],z['c'],z['d'])
        oldrows=sorted(key(z) for z in old['last_unit_prefixes']);newrows=sorted((z['D'],z['a_m'],z['a_n'],z['b_m'],z['b_n'],z['c'],z['d']) for z in unit)
        assert len(set(newrows))==194 and oldrows==newrows
    out={'status':'passed_recomputed_all_396738_lattice_pairs_and_fourth_row_exclusion','scope':'Only degree7 real quadratic P with nonrational positions {1,2}, or reversal {5,6}; arbitrary finite integer cofactors. General conjecture open.','low_counts':report['counts'],'cd_counts':dict(counts),'unit_cases':194,'all_C4':2,'survivors':0,'elapsed_seconds':time.monotonic()-begin,'python':sys.version,'manifest_sha256':sha((r/'manifest.json').read_bytes()),'verifier_sha256':sha(Path(__file__).read_bytes())}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
