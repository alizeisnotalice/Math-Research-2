"""Rebuild the .83743 whole-leaf seed by auditing the standard and signed score indices."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent


def sha(data):
    return sha256(data).hexdigest()


def main():
    if not __debug__:
        raise SystemExit('Assertions are required.')
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('computation', 'whole_index_083740', 'whole_index_083742', 'work'):
        parser.add_argument(key, type=Path)
    parser.add_argument('--replay-standard', nargs=3, type=Path, metavar=('CERT', 'TABLE32', 'CLANG_CACHE'))
    parser.add_argument('--workers', type=int, default=4)
    a = parser.parse_args()
    assert 1 <= a.workers <= 12
    main, w = a.computation.resolve(), a.work.resolve()
    w.mkdir(exist_ok=True, parents=True)
    stage = w/'standard-stage'
    pins = json.loads((HERE/'tomar-cross-index-083743-provenance.json').read_bytes())['inputs']

    def run(script, args, outfile=None):
        command = [sys.executable, '-B', str(script), *map(str, args)]
        if outfile is None:
            subprocess.run(command, check=True)
        else:
            with outfile.open('wb') as output:
                subprocess.run(command, check=True, stdout=output)

    bundle = (HERE/'tomar-standard-083741-reuse.zip').read_bytes()
    assert sha(bundle) == pins['standard41_bundle']
    with ZipFile(BytesIO(bundle)) as package:
        ledger = package.read('stage-ledger.json')
        audit = package.read('standard41-audit.json')
        expected_standard = package.read('standard41-index.json')
        assert sha(ledger) == pins['standard41_ledger']
        assert sha(audit) == pins['standard41_audit']
        assert sha(expected_standard) == pins['standard41_complete']
        stage.mkdir(exist_ok=True)
        (stage/'stage-ledger.json').write_bytes(ledger)
        (stage/'retained.json').write_text(json.dumps(json.loads(audit)['retained']))
        parts = {p for p,h in json.loads(ledger)['audits'][0]['work_files']}
        assert len(parts) == 1024
        expected = {f'stage{k}/part{p:04}.{ext}' for p in parts for k in range(1,5)
                    for ext in ('json','npz')}
        seen = set()
        for j in (0,1):
            raw = package.read(f'tomar-scored-083741-{j}.zip')
            assert sha(raw) == pins[f'standard41_archive_{j}']
            with ZipFile(BytesIO(raw)) as shard:
                for name in shard.namelist():
                    assert name in expected and name not in seen
                    seen.add(name)
                    dest = stage/name
                    dest.parent.mkdir(exist_ok=True)
                    dest.write_bytes(shard.read(name))
        assert seen == expected

    base40 = main/'tomar-multi-target-after-083740.json'
    base42 = main/'tomar-multi-target-after-083742.json'
    radial = main/'five-body-radial-affine.json'
    assert sha(base40.read_bytes()) == pins['standard40_complete']
    assert sha(base42.read_bytes()) == pins['signed42_complete']
    for code, target in [('41','.83741'), ('42','.83742')]:
        dest = w/('base'+code+'.json')
        run(main/'prepare_tomar_multi_target_layer.py',
            [base40, a.whole_index_083740, radial, w/('base'+code), '--target', target,
             '--radial-target', '.83742'], dest)
        expected_seed = main/('tomar-083741-seed.json' if code=='41' else 'tomar-083742-future-seed.json')
        assert dest.read_bytes() == expected_seed.read_bytes()
    frozen = json.loads(ledger)
    if a.replay_standard:
        tree, table, cache = a.replay_standard
        fresh = w/'standard-geometric-replay'
        run(main/'replay_tomar_083731_scored.py',
            [main, tree, table, w/'base41', w/'base41.json', fresh, cache, '--workers', a.workers])
        actual = json.loads((fresh/'stage-ledger.json').read_bytes())
        assert {k:v for k,v in frozen.items() if k != 'binary_sha'} == {
            k:v for k,v in actual.items() if k != 'binary_sha'}
        for name in expected:
            assert (fresh/name).read_bytes() == (stage/name).read_bytes()
    run(main/'audit_tomar_083731_score_layers.py',
        [stage,w/'base41',w/'base41.json',w/'base42',w/'base42.json'],w/'standard-audit.json')
    actual_audit = json.loads((w/'standard-audit.json').read_bytes())
    old_audit = json.loads(audit)
    for key in actual_audit:
        assert actual_audit[key] == old_audit[key]
    (w/'canonical-standard-audit.json').write_bytes(audit)
    run(main/'merge_tomar_multi_target_layer.py',
        [base40,a.whole_index_083740,stage,w/'base41.json',w/'canonical-standard-audit.json',
         w/'standard-index'],w/'standard-index.json')
    assert (w/'standard-index.json').read_bytes() == expected_standard
    combine = HERE/'combine_tomar_whole_leaf_indices.py'
    assert sha(combine.read_bytes()) == pins['pointwise_max_source']
    run(combine,[main,base42,a.whole_index_083742,w/'standard-index.json',w/'standard-index',
                 main/'tomar-metal-integer-floors-083742.json',main/'tomar-metal-integer-floors-083742.npz',
                 w/'combined'],w/'combined.json')
    assert sha((w/'combined.json').read_bytes()) == pins['combined_complete']
    for code, target in [('43','.83743'),('44','.83744')]:
        dest = w/('seed'+code+'.json')
        run(main/'prepare_tomar_multi_target_layer.py',
            [w/'combined.json',w/'combined',radial,w/('seed'+code),'--target',target,
             '--radial-target','.83744'],dest)
        assert sha(dest.read_bytes()) == pins['seed'+code]
    print(json.dumps(dict(complete=True,seed43=json.loads((w/'seed43.json').read_bytes())['remaining'],
                          seed44=json.loads((w/'seed44.json').read_bytes())['remaining'])))


if __name__ == '__main__':
    main()
