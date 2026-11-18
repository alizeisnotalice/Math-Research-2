"""Generate and verify the 36 original ordinary batches, including all 18 depth stops."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('computation', 'tree', 'table32', 'table64', 'work', 'cache'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert __debug__ and 1 <= args.workers <= 32
    args.work.mkdir(parents=True, exist_ok=False)
    args.cache.mkdir(parents=True, exist_ok=True)
    package = Path(__file__).resolve().parent
    with tarfile.open(package/'volumes/native45-receipts.tar.gz') as ar:
        manifest = {name: (size, digest) for name, size, digest in json.load(ar.extractfile('MANIFEST.json'))}

        def read(name):
            raw = ar.extractfile(name).read()
            assert (len(raw), sha256(raw).hexdigest()) == manifest[name]
            return raw

        code = args.work/'produce_p30_083745.py'
        for path in (code, code.with_name('finite-point-lift-083744.json')):
            path.write_bytes(read('inputs/ordinary/'+path.name))
        batch_receipts = sorted(name for name in manifest if name.startswith('ordinary/') and name.endswith('/result.json'))
        assert len(batch_receipts) == 36
        seen = set()
        failed = set()
        total_nodes = 0
        for receipt_name in batch_receipts:
            prefix = receipt_name.removesuffix('result.json')
            group = prefix.split('/')[1]
            selection = args.work/('group-'+group.rsplit('-', 1)[1]+'.json')
            raw = read('inputs/ordinary/'+selection.name)
            selection.write_bytes(raw)
            keys = {tuple(pair) for pair in json.loads(raw)['leaves']}
            assert not (keys & seen)
            seen.update(keys)
            dest = args.work/group
            cmd = [sys.executable, '-I', '-B', code, args.computation, args.tree,
                   args.table32, args.table64, dest, args.cache, selection, '--workers', str(args.workers)]
            result = json.loads(subprocess.check_output(cmd))
            expected = json.loads(read(receipt_name))
            assert result['trees'] == expected['trees'] == len(keys)
            assert result['selection_sha256'] == expected['selection_sha256'] == sha256(raw).hexdigest()
            failures = set(map(tuple, result['failed']))
            assert failures == set(map(tuple, expected['failed']))
            failed.update(failures)
            for part, pos in keys:
                filename = f'part{part:04}-{pos}.json'
                row = json.loads((dest/filename).read_bytes())
                assert row == json.loads(read(prefix+'rows/'+filename))
                total_nodes += row['stat'][1]
            print(f'{group} verified, {len(seen)} leaves checked', flush=True)
    assert (len(seen), len(failed), total_nodes) == (3071, 18, 1161865226)
    record = dict(parts=36, ordinary_attempts=len(seen), success=len(seen)-len(failed),
                  depth_stops=sorted(failed), final_attempt_nodes=total_nodes)
    (args.work/'summary.json').write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    main()
