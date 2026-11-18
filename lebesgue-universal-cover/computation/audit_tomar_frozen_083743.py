"""Reconstruct the authenticated .83743 seed and audit all four final packet layers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('computation','index40','index42','work'):
        parser.add_argument(key,type=Path)
    a=parser.parse_args()
    computation=a.computation.resolve()
    work=a.work.resolve()
    subprocess.run([sys.executable,'-I','-B',HERE/'replay_tomar_cross_index_083743.py',
                    computation,a.index40,a.index42,work/'index'],check=True)
    report_file=HERE/'tomar-083743-signed-screen-audit.json'
    blob=report_file.read_bytes()
    assert sha256(blob).hexdigest()=='890f95851261047bf8364476626129cad6fe82b832994e3ffd8047f318ae75eb'
    report=json.loads(blob)
    for name,h in report['files'].items():
        assert sha256((HERE/name).read_bytes()).hexdigest()==h
    stage=work/'stage'
    stage.mkdir(exist_ok=True)
    ledger=(HERE/'tomar-083743-stage-ledger.json').read_bytes()
    (stage/'stage-ledger.json').write_bytes(ledger)
    (stage/'retained.json').write_bytes((HERE/'tomar-083743-retained.json').read_bytes())
    actual=json.loads(ledger)
    members={f'stage{k}/part{p:04}.{suffix}'
             for k,layer in enumerate(actual['audits'],1) for p,h in layer['work_files']
             for suffix in ('npz','json')}
    seen=set()
    for archive in report['archives']:
        file=HERE/archive['bundle']
        assert sha256(file.read_bytes()).hexdigest()==archive['sha256']
        with ZipFile(file) as z:
            for name in z.namelist():
                assert name in members and name not in seen
                seen.add(name)
                out=stage/name
                out.parent.mkdir(exist_ok=True)
                out.write_bytes(z.read(name))
    assert len(seen)==8192 and seen==members
    with (work/'audit.json').open('wb') as output:
        subprocess.run([sys.executable,'-B',computation/'audit_tomar_083731_score_layers.py',
                        stage,work/'index/seed43',work/'index/seed43.json',work/'index/seed44',
                        work/'index/seed44.json'],check=True,stdout=output)
    assert (work/'audit.json').read_bytes()==(HERE/'tomar-083743-score-layers-audit.json').read_bytes()
    print(json.dumps(dict(complete=True,parts=1024,files=len(seen),
                          remaining=report['remaining'],stage_calls=sum(s['calls']for s in report['stages']))))


if __name__=='__main__':
    main()
