"""Rebuild the exact tables and seed from the archived geometry and public project sources."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('computation','tree','work'):p.add_argument(name,type=Path)
    a=p.parse_args();assert __debug__
    comp=a.computation.resolve();tree=a.tree.resolve();work=a.work.resolve();work.mkdir(parents=True,exist_ok=False)
    package=Path(__file__).resolve().parent
    sha=lambda path:sha256(path.read_bytes()).hexdigest()
    index=work/'index';index.mkdir();report=comp/'tomar-multi-target-after-083743.json'
    assert sha(tree/'SHA256SUMS')=='56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410'
    assert sha(report)=='6f4df8e2a40c0a6a0f88ad0b7688b8a5ff7459dff1d02283aeb9867f88027459'
    files={f'part{r[0]:04}.npz' for r in json.loads(report.read_bytes())['rows']};assert len(files)==1024
    with ZipFile(package/'volumes/native45-prior-index.zip') as z:
        assert set(z.namelist())==files and len(z.namelist())==1024
        for f in sorted(files): (index/f).write_bytes(z.read(f))
    sys.path.insert(0,str(comp))
    from prepare_tomar_multi_target_layer import prepare
    seed=prepare(argparse.Namespace(index_report=report,index=index,radial_certificate=comp/'five-body-radial-affine.json',emit=work/'seeds',target=F('.83745'),radial_target=F('.83745')))
    raw=(json.dumps(seed,separators=(',',':'))+'\n').encode();assert raw==(comp/'tomar-083745-future-seed.json').read_bytes();(work/'seed.json').write_bytes(raw)
    def build(n):
        table=work/f'table{n}.npz';subprocess.run([sys.executable,'-I','-B',comp/'tomar_fan_table.py',tree,table,'--subdiv',str(n)],check=True)
        return sha(table)
    with ThreadPoolExecutor(2) as ex: hashes=dict(zip((32,64),ex.map(build,(32,64))))
    assert hashes=={32:'fe5b816060d713642dc44e3f8981f4ae5d1b2b1c1e271b7dff75f031dd75d87c',64:'0913f79e54bfcb804c8c4a4c0bbd7962d42d461591175ab037739b3562c6016d'}
    print(json.dumps({'seed':sha(work/'seed.json'),'parts':len(files),'tables':hashes}))

if __name__=='__main__':main()
