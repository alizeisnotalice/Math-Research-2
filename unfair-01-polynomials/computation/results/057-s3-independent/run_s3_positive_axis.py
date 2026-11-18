"""Verify compressed interval shards and evaluate exact positivity counts."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

BINARY = Path(__file__).resolve().parent / 's3_positive_axis'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run_shard(root, entry, stride):
    path = (root / entry['path']).resolve()
    assert path.is_relative_to(root) and not path.is_symlink()
    compressed = path.read_bytes()
    assert len(compressed) == entry['bytes'] and sha(compressed) == entry['sha256']
    raw = gzip.decompress(compressed)
    assert len(raw) == entry['raw_bytes'] and sha(raw) == entry['raw_sha256']
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == entry['groups']
    assert len({tuple(row[5:7]) for row in rows}) == len(rows)
    for row in rows:
        assert len(row) == 14 and all(type(value) is int for value in row)
        assert row[:5] == entry['core']
        assert math.isqrt(row[4])**2 != row[4] and (row[2]**2 - row[4]) % 4 == 0
        assert -120 <= row[5] <= 120 and -92 <= row[6] <= 92
    plain = '\n'.join(' '.join(map(str, row)) for row in rows) + '\n'
    stats = json.loads(subprocess.check_output([str(BINARY), str(stride)], input=plain, text=True))
    assert stats['groups'] == entry['groups'] and stats['mothers'] == entry['mothers']
    return {'core': entry['core'], 'raw_sha256': entry['raw_sha256'], **stats}


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    p = ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--stride', type=int, default=257)
    p.add_argument('--expect-mothers', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert not a.out.exists() and 1 <= a.workers <= 13 and a.stride >= 1
    manifest = a.manifest.read_bytes()
    entries = [json.loads(line) for line in manifest.splitlines()]
    assert len({tuple(x['core']) for x in entries}) == len(entries)
    nonempty = [x for x in entries if x['path'] is not None]
    assert len({x['path'] for x in nonempty}) == len(nonempty)
    assert all(x['groups']==x['mothers']==x['bytes']==x['raw_bytes']==0 for x in entries if x['path'] is None)
    assert sum(x['mothers'] for x in entries) == a.expect_mothers
    root = a.manifest.parent.resolve()
    counters = Counter()
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        rows = list(pool.map(lambda x: run_shard(root, x, a.stride), nonempty))
    for row in rows:
        counters.update({k: v for k,v in row.items() if type(v) is int})
    assert counters['mothers'] == a.expect_mothers and counters['axis_pass'] + counters['axis_reject'] == counters['mothers']
    out = {'status':'verified_integer_axis_plane_counts', 'manifest_sha256':sha(manifest),
           'source_sha256':sha(BINARY.with_suffix('.cpp').read_bytes()), 'binary_sha256':sha(BINARY.read_bytes()),
           'cores':len(entries), 'nonempty_cores':len(nonempty), 'workers':a.workers, 'stride':a.stride, 'counts':dict(counters),
           'elapsed_seconds':time.monotonic()-started, 'core_records':rows,
           'scope':'Necessary positivity on nine rational points only. Complete coverage refers to the bound input ledger.'}
    with a.out.open('x') as f:
        f.write(json.dumps(out, separators=(',',':'))+'\n')
    print(json.dumps({k:out[k] for k in ('status','cores','counts','elapsed_seconds')}))


if __name__ == '__main__':
    main()
