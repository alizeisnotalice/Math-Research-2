"""Recheck each root and leaf at an increased exact decimal threshold."""
import argparse, concurrent.futures, hashlib, json, os, subprocess, sys, tempfile
from pathlib import Path

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','directory','scratch','receipt'):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    original=json.loads(a.receipt.read_text())
    assert original['complete'] and original['parts']==list(range(16))
    manifest=(a.source/'certificates/SHA256SUMS').read_text()
    expected='2026ebe7e5339b89542e0429e5ad704b1a956278e3ae455b6d667b4c11892846'
    archive=a.source/'certificates/cert5d_08347.tar.xz'
    assert digest(archive)==expected==original['archive_sha256'] and expected in manifest
    import tarfile
    with tarfile.open(archive) as tar:
        embedded=next(x for x in tar if x.name.endswith('/SHA256SUMS'))
        records=tar.extractfile(embedded).read().decode().splitlines()
    checker=a.source/'checker/lbcheck.py'
    assert digest(checker)==original['checker_sha256']['lbcheck.py']
    target='0.834700003'
    with tempfile.TemporaryDirectory(dir=a.scratch) as temp:
        root=Path(temp)
        for record in records:
            sha,name=record.split();p=a.directory/name
            assert digest(p)==sha
            if p.suffix=='.json':
                header=json.loads(p.read_text())
                assert str(header['tau'])=='0.8347' and header['ns']==[3,5]
                header['tau']=float(target); assert repr(header['tau'])==target
                (root/p.name).write_text(json.dumps(header))
            else:
                (root/p.name).symlink_to(p.resolve())
        def check(i):
            env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
            raw=subprocess.check_output([sys.executable,'-I','-B',str(checker),str(root),str(i)],env=env,text=True)
            d=json.loads(raw);assert d['part']==i and d['ok'] and d['fails']==0 and d['leaves']==original['trees'][i]['leaves']
            return d
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            rows=list(pool.map(check,range(16)))
    print(json.dumps(dict(target=target,source_archive=expected,checker=digest(checker),
        leaves=sum(z['leaves'] for z in rows),parts=rows),indent=2))
if __name__=='__main__':main()
