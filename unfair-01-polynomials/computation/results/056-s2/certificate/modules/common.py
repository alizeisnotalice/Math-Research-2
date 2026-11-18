"""Deterministic certificate framing and safe readers."""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE_CODE = {
    'stage4': 's2_degree8_stage4.py',
    'stage58': 's2_degree8_stage5to8.py',
    'root': 's2_degree8_root_screen.py',
    'trees': 's2_degree8_dual_trees.py',
}
CORE_STAGES = ('stage4', 'stage58', 'root')
EXPECTED_CORE_HASH = 'f3668066ae89d45115788b860d8093d81c6b90e6cae9116ddaa7ea0a3db8b45b'
EXPECTED_TREE_HASH = 'bcd1f36bdeb1499edf9528d50d72227adefc68569e7a7f9846f3cb4c6b1e5421'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def guard_assertions():
    if sys.flags.optimize:
        raise SystemExit('refuse optimized Python: the production modules require assertions')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def canonical(row):
    return json.dumps(row, sort_keys=True, separators=(',', ':')).encode() + b'\n'


def payload_hash(rows):
    h = hashlib.sha256()
    for row in rows:
        h.update(canonical(row))
    return h.hexdigest()


def read_json(path):
    return json.loads(path.read_bytes())


def input_hashes(data_path):
    gzip_path = data_path.with_suffix('.jsonl.gz')
    meta = read_json(data_path.with_suffix('.json'))
    raw = gzip.decompress(gzip_path.read_bytes())
    rows = [json.loads(line) for line in raw.splitlines()]
    require(meta['status'] == 'complete' and meta['mothers'] == len(rows) == 1583, 'tree input count')
    require(digest(raw) == meta['rows_payload_sha256'] == EXPECTED_TREE_HASH, 'tree input payload binding')
    require(meta['gzip_sha256'] == sha(gzip_path), 'tree input gzip binding')
    require(rows == sorted(rows) and len({tuple(row) for row in rows}) == len(rows), 'tree input order/uniqueness')
    return rows, meta


def replay_files(out, stage):
    parent = out / stage
    fileglob = 'mother-*.json' if stage == 'trees' else '*.json'
    base = parent if stage == 'trees' else parent / 'cores'
    for path in sorted(base.glob(fileglob)):
        yield path.name, path.read_bytes()


def capsule_files(index, stage):
    previous = None
    for part in index['parts'][stage]:
        path = ROOT / 'data' / part['path']
        count = 0
        first = None
        with tarfile.open(path, mode='r|gz') as archive:
            for item in archive:
                name = item.name
                bare = Path(name).name == name and name.endswith('.json')
                require(item.isreg() and bare and (previous is None or previous < name), f'{stage}: unsafe/unordered {name}')
                if first is None:
                    first = name
                count += 1
                previous = name
                yield name, archive.extractfile(item).read()
        require(count == part['count'] and first == part['first'] and previous == part['last'], f'{stage}: TOC mismatch')
