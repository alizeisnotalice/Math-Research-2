#!/usr/bin/env python3
"""Check byte bindings and coverage of all corrected S={2} replay records."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True
from common import (
    CORE_STAGES, EXPECTED_CORE_HASH, ROOT, STAGE_CODE, canonical, capsule_files, digest,
    guard_assertions, input_hashes, payload_hash, read_json, replay_files, require, sha,
)

EMPTY_STATE_HASH = digest(b'')
CORE_COUNT = 11760
STAGE4_COUNT = 219
MOTHER_COUNT = 461683
HASH_FORMAT = re.compile(r'[0-9a-f]{64}\Z')
POINTS = [[1, 1], [1, 2], [2, 1], [1, 4], [4, 1], [3, 4], [4, 3], [6, 5], [5, 6]]


def verify_manifest(root):
    manifest = read_json(root / 'manifest.json')
    require(manifest['format'] == 2, 'manifest format')
    wanted = set()
    for row in manifest['files']:
        rel = Path(row['path'])
        require(not rel.is_absolute() and '..' not in rel.parts and rel not in wanted, f'manifest path {rel}')
        path = root / rel
        require(path.is_file() and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], f'manifest bytes {rel}')
        wanted.add(rel)
    listed = list(root.rglob('*'))
    require(not any(path.is_symlink() for path in listed), 'symlinks in package')
    actual = {path.relative_to(root) for path in listed if path.is_file() and path != root / 'manifest.json'}
    require(actual == wanted, f'manifest inventory mismatch: {actual ^ wanted}')
    return sha(root / 'manifest.json')


def core_records(files, stage, summary, module_hash):
    records_hash, raw_hash = hashlib.sha256(), hashlib.sha256()
    cores = {}
    counts = Counter()
    size = 0
    previous = None
    for name, raw in files:
        node = json.loads(raw)
        bind = node['binding']
        core = tuple(bind['core'])
        require(node['status'] == 'complete' and bind['code_sha256'] == module_hash, f'{stage}: stale or incomplete {name}')
        require((previous is None or previous < name) and len(core) == 4 and core not in cores, f'{stage}: duplicate or unordered core {name}')
        previous = name
        L, T, N, delta = core
        require(name == f'L{L}_T{T}_N{N}.json' and delta == T*T-4*N, f'{stage}: core identity {name}')
        rows = node['records']
        stats = node['summary']
        key = {'stage4': 'row4_states', 'stage58': 'row8_states', 'root': 'mothers'}[stage]
        require(stats[key] == len(rows), f'{stage}: record count {name}')
        if stage == 'root':
            require(tuple(stats['core']) == core and bind['points'] == POINTS, f'{stage}: metadata {name}')
            actual = Counter()
            for row in rows:
                require(tuple(row['row'][:4]) == core and len(row['row']) == 9
                        and len(row['H']) == 17 and row['H'][0] == row['H'][-1] == 1, f'{stage}: row shape {name}')
                status = 'passed' if row['code'] == 255 else 'rejected'
                require(row['status'] == status, f'{stage}: row code {name}')
                actual[status] += 1
                actual[f"code_{row['code']}"] += 1
            require(actual == stats['counts'], f'{stage}: row status counts {name}')
        else:
            require((stats['L'], stats['T'], stats['N'], stats['Delta']) == core, f'{stage}: metadata {name}')
            width = 12 if stage == 'stage4' else 21
            require(all(len(row) == width and tuple(row[:4]) == core for row in rows), f'{stage}: unbound row {name}')
        require(len(rows) == len({canonical(row) for row in rows}), f'{stage}: duplicate record {name}')
        counts.update(stats['counts'])
        size += len(rows)
        raw_hash.update(name.encode() + b'\0' + raw)
        records_hash.update(canonical([core, stats['counts'], rows]))
        cores[core] = (digest(raw), bind, rows)
    expected_sizes = {'stage4': STAGE4_COUNT, 'stage58': MOTHER_COUNT, 'root': MOTHER_COUNT}
    require(len(cores) == summary['cores'] == summary['completed'] == CORE_COUNT and summary['status'] == 'complete'
            and size == summary[key] == expected_sizes[stage]
            and dict(counts) == summary['counts'] and summary['code_sha256'] == module_hash, f'{stage}: aggregate summary')
    require(payload_hash(sorted(cores)) == EXPECTED_CORE_HASH, f'{stage}: full key universe')
    return cores, raw_hash.hexdigest(), records_hash.hexdigest()


def audit_trees(files, rows, source, summary, module_hash):
    tree_hash = hashlib.sha256()
    depths, live_maximums = Counter(), Counter()
    seen = set()
    count = 0
    settings = {'depth': 128, 'cap': 25000, 'seconds_per_tree': 20.0}
    require(summary['status'] == 'complete' and summary['code_sha256'] == module_hash
            and summary['input_payload_sha256'] == source['rows_payload_sha256']
            and summary['depth'] == settings['depth'] and summary['cap'] == settings['cap']
            and summary['tree_seconds'] == settings['seconds_per_tree'], 'tree summary settings')
    for name, raw in files:
        result = json.loads(raw)
        row = result['binding']['row']
        key = tuple(row)
        file_id = digest(canonical(row)[:-1])
        check = {'version': 1, 'row': row, 'input_payload_sha256': source['rows_payload_sha256'],
                 'code_sha256': module_hash, **settings}
        require(name == f'mother-{file_id}.json' and key not in seen and result['binding'] == check
                and result['status'] == 'complete', f'tree identity {name}')
        seen.add(key)
        ab, hgq = result['ab'], result['hgq']
        require(ab == hgq and ab['status'] == 'dead_all_prefixes', f'trees differ or have unknown state {name}')
        trace = ab['trace']
        require(1 <= len(trace) == ab['depth_reached'] <= settings['depth'], f'tree depth {name}')
        for i, step in enumerate(trace, 1):
            require(step['k'] == i and 0 <= step['nonbinary_C'] and HASH_FORMAT.fullmatch(step['state_sha256'])
                    and (step['live'] > 0 if i < len(trace) else step['live'] == 0)
                    and (i < len(trace) or step['state_sha256'] == EMPTY_STATE_HASH), f'tree layer {name} {i}')
        count += 1
        depths[len(trace)] += 1
        live_maximums[max(step['live'] for step in trace)] += 1
        tree_hash.update(name.encode() + b'\0' + raw)
    require(seen == {tuple(row) for row in rows} and count == len(rows) == summary['mothers'] == summary['completed']
            and summary['counts'] == {'tree_dead_all_prefixes': count}, 'incomplete dual-tree mother coverage')
    return {'dead': count, 'unknown': 0, 'death_depths': dict(sorted(depths.items())),
            'max_live_distribution': dict(sorted(live_maximums.items())),
            'state_layers': sum(k*v for k, v in depths.items()), 'records_sha256': tree_hash.hexdigest()}


def audit_run(output=None):
    index = read_json(ROOT / 'data/index.json') if output is None else None
    modules = {stage: sha(ROOT / 'modules' / filename) for stage, filename in STAGE_CODE.items()}
    if output is None:
        read_stage = lambda stage: capsule_files(index, stage)
        summaries = {stage: ROOT / 'frozen' / f'{stage}-summary.json' for stage in STAGE_CODE}
        input_prefix = ROOT / 'frozen/tree-input'
    else:
        read_stage = lambda stage: replay_files(output, stage)
        summaries = {stage: output / stage / 'summary.json' for stage in STAGE_CODE}
        input_prefix = output / 'tree-input'
    summary = {stage: read_json(path) for stage, path in summaries.items()}
    hashes = {stage: sha(path) for stage, path in summaries.items()}
    tables, raw_hashes, semantic_hashes = {}, {}, {}
    for stage in CORE_STAGES:
        tables[stage], raw_hashes[stage], semantic_hashes[stage] = core_records(read_stage(stage), stage, summary[stage], modules[stage])
    stage4, stage58, root = (tables[stage] for stage in CORE_STAGES)
    join_hash = digest((hashes['stage4'] + raw_hashes['stage4']).encode())
    require(set(stage4) == set(stage58) == set(root) and summary['stage58']['stage4_input_sha256'] == join_hash
            and summary['root']['stage58_input_sha256'] == raw_hashes['stage58'], 'stages summary binding')
    code255 = []
    for core, (_, binding, records) in root.items():
        source_sha58, binding58, source_rows58 = stage58[core]
        source_sha4, _, source_rows4 = stage4[core]
        require(binding58['stage4_source_sha256'] == source_sha4 and binding58['stage4_input_sha256'] == join_hash
                and binding['stage58_source_sha256'] == source_sha58
                and binding['stage58_input_sha256'] == raw_hashes['stage58'], f'core dependency binding {core}')
        ledger4 = {tuple(row) for row in source_rows4}
        require(all(tuple(row[:6] + row[9:15]) in ledger4 for row in source_rows58), f'stage58: unknown parent {core}')
        rows58 = {tuple(row[:9]) for row in source_rows58}
        require(len(rows58) == len(source_rows58) == len(records)
                and [tuple(record['row']) for record in records] == sorted(rows58), f'root: incomplete mother list {core}')
        code255.extend(record['row'] for record in records if record['code'] == 255)
    rows, input_meta = input_hashes(input_prefix)
    require(input_meta['source_summary_sha256'] == hashes['root'] and input_meta['source_core_ledger_sha256'] == raw_hashes['root']
            and input_meta['source_code_sha256'] == modules['root'] and rows == sorted(code255), 'tree input not all code-255 rows')
    trees = audit_trees(read_stage('trees'), rows, input_meta, summary['trees'], modules['trees'])
    semantic = {'stage_records_sha256': semantic_hashes, 'tree_input_sha256': input_meta['rows_payload_sha256'],
                'stage_counts': {stage: summary[stage]['counts'] for stage in CORE_STAGES}, 'trees': trees}
    return {'status': 'complete', 'semantics': semantic, 'raw_core_ledger_sha256': raw_hashes, 'summary_sha256': hashes}


def main():
    guard_assertions()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path, help='also check replay output and compare complete deterministic commitments')
    options = parser.parse_args()
    manifest_hash = verify_manifest(ROOT)
    original = audit_run()
    report = {'status': 'passed_archive_binding', 'manifest_sha256': manifest_hash, 'original': original}
    if options.replay is not None:
        replay = audit_run(options.replay)
        require(replay['semantics'] == original['semantics'], 'full replay differs from sealed records')
        report.update(status='passed_archive_and_replay_binding', replay=replay)
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
