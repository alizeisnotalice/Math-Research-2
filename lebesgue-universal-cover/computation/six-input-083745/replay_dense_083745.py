"""Rebuild the 24 complete mesh64-dense trees, or verify an existing independent run."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tarfile


def inputs(package):
    jobs = []
    for filename, prefix, result_prefix in (
        ('native45-receipts.tar.gz', 'inputs/dense/', 'dense/native45-dense-'),
        ('native45-pure-a.tar.gz', 'dense/', 'dense/proof-'),
    ):
        with tarfile.open(package / filename) as archive:
            manifest = {name: (size, digest) for name, size, digest in json.load(archive.extractfile('MANIFEST.json'))}

            def checked(name):
                raw = archive.extractfile(name).read()
                assert (len(raw), sha256(raw).hexdigest()) == manifest[name]
                return raw

            for name in manifest:
                if not name.startswith(prefix) or not name.endswith('.json') or '/' in name[len(prefix):]:
                    continue
                selection = checked(name)
                data = json.loads(selection)
                key = f"{data['part']}-{data['pos']}"
                assert key == Path(name).stem and data['paths'] == [''] and data['terminal_guard'] == 'mesh64-dense'
                receipt = checked(result_prefix + key + '/result.json')
                expected = json.loads(receipt)
                assert sha256(selection).hexdigest() == expected['selection_sha256']
                assert expected['all_closed'] and expected['complete'] and expected['threshold_integer'] == 120689264254125501
                jobs.append((key, selection, expected))
    assert len(jobs) == len({job[0] for job in jobs}) == 24
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('computation', 'tree', 'table64', 'runs', 'cache', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--check-existing', action='store_true', help='Compare existing result.json files without running the area kernels.')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    assert __debug__ and args.workers >= 1
    jobs = inputs(Path(__file__).resolve().parent / 'volumes')
    if not args.check_existing:
        args.runs.mkdir(parents=True, exist_ok=False)
        args.cache.mkdir(parents=True, exist_ok=True)

    def run(job):
        key, selection, expected = job
        run_dir = args.runs / key
        if not args.check_existing:
            selection_path = args.runs / (key + '.json')
            selection_path.write_bytes(selection)
            cmd = [sys.executable, '-I', '-B', args.computation / 'replay_scaled_p30_high_leaf.py',
                   args.computation, args.tree, args.table64, run_dir, args.cache, selection_path]
            with (args.runs / (key + '.log')).open('w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        raw = (run_dir / 'result.json').read_bytes()
        actual = json.loads(raw)
        assert actual == expected, key
        return dict(key=key, nodes=actual['rows'][0]['stat'][1], receipt_sha256=sha256(raw).hexdigest())

    jobs.sort(key=lambda job: -job[2]['rows'][0]['stat'][1])
    with ThreadPoolExecutor(args.workers) as workers:
        rows = sorted(workers.map(run, jobs), key=lambda row: row['key'])
    total = sum(row['nodes'] for row in rows)
    assert total == 32803666
    report = dict(count=len(rows), nodes=total, all_full_field_equal_and_closed=True, cases=rows)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(f"{len(rows)} / {len(jobs)} trees closed, {total} nodes")


if __name__ == '__main__':
    main()
