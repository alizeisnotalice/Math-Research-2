#!/usr/bin/env python3
"""Verify frozen packet and recovered-result bytes; this is not a proof checker."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(path, row):
    assert path.is_file(), path
    actual = sha(path)
    assert actual == row['sha256'], (path, actual, row['sha256'])
    if 'bytes' in row:
        assert path.stat().st_size == row['bytes'], path
    return actual


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spec', type=Path, required=True)
    ap.add_argument('--packet', type=Path, action='append', required=True)
    ap.add_argument('--report', type=Path, required=True)
    ap.add_argument('--review-root', type=Path, required=True, help='topic directory')
    ap.add_argument('--archive-root', type=Path, help='repository root; use archived destinations rather than producer paths')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    spec = json.loads(args.spec.read_text())
    assert len(args.packet) == len(spec['packets'])
    packet_results, resolved = [], {}
    for packet in args.packet:
        expected = spec['packets'][packet.name]
        check(packet, expected)
        data = json.loads(packet.read_text())
        topic = packet.parents[2]
        rows = []
        for row in data['files'] + data.get('existing_dependencies', []):
            dest = row.get('destination', 'planar-centered-square-weak11/' + row.get('path', ''))
            path = args.archive_root / dest if args.archive_root else (Path(row['source']) if 'source' in row else topic / row['path'])
            check(path, row)
            resolved[dest] = path
            rows.append(dict(destination=dest, bytes=path.stat().st_size, sha256=row['sha256']))
        inherited = data.get('inherited_source', {})
        inherited_rows = inherited.get('paths', [])
        if 'inherited_fixture' in data:
            inherited_rows = inherited_rows + [data['inherited_fixture']]
        for row in inherited_rows:
            dest = 'planar-centered-square-weak11/' + row['path']
            path = args.archive_root / dest if args.archive_root else Path(spec['original_binding_paths'][dest])
            check(path, row)
            rows.append(dict(destination=dest, bytes=path.stat().st_size, sha256=row['sha256'], role='inherited source binding only'))
        packet_results.append(dict(manifest_destination=expected['destination'], manifest_sha256=expected['sha256'], checked_files=rows))
    check(args.report, spec['report'])
    recovered = []
    for item in spec['recovered_results']:
        path = args.review_root / item['path']
        check(path, item)
        data = json.loads(path.read_text())
        assert data['source_sha256'] == sha(args.review_root / item['source'])
        for field, relative in item.get('embedded_bindings', {}).items():
            target = args.archive_root / relative if args.archive_root else Path(spec['original_binding_paths'][relative])
            assert data[field] == sha(target), (field, target)
        recovered.append(dict(path=item['path'], sha256=item['sha256'], source_sha256=data['source_sha256'], python=data['python'], elapsed_seconds=data['elapsed_seconds'], execution_receipt='not persisted; original exit-zero claim only in contemporaneous execution report, not re-observed'))
    # The saved production replay pins the compute packet and all its payloads/dependencies.
    leaf = json.loads((args.review_root / 'computation/results/round070-independent/saved-cells.json').read_text())
    compute = next(p for p in packet_results if p['manifest_destination'].endswith('round070-compute-integration-manifest.json'))
    assert leaf['manifest_sha256'] == compute['manifest_sha256']
    assert leaf['pins'] == {r['destination']: r['sha256'] for r in compute['checked_files']}
    base = 'planar-centered-square-weak11/computation/'
    prefix_receipt = json.loads(resolved[base+'results/round070-root-prefix-grid/grid-prefix-receipt.json'].read_text())
    assert prefix_receipt['exit_status'] == 0 and prefix_receipt['workers_requested'] == 12
    assert prefix_receipt['distinct_workers_with_results'] == 8
    assert len(prefix_receipt['owned_pids_absent']) == 9 and all(prefix_receipt['owned_pids_absent'].values())
    assert prefix_receipt['result_sha256'] == sha(resolved[base+'results/round070-root-prefix-grid/r070-grid-prefix-result.json'])
    wl_receipt = json.loads(resolved[base+'results/round070-root-symbolic/wolfram-receipt.json'].read_text())
    wl_result = json.loads(resolved[base+'results/round070-root-symbolic/r070-captured-mass-symbolic.json'].read_text())
    assert wl_receipt['exit_status'] == 0 and wl_receipt['owned_pid_absent']
    assert wl_receipt['source_sha256'] == sha(resolved[base+'round070/r070-captured-mass-symbolic.wl'])
    assert len(wl_result['checks']) == 14 and all(wl_result['checks'].values())
    finite_result = args.review_root / 'computation/results/round070-independent/finite-prefix-area.json'
    finite_receipt = args.review_root / 'computation/results/round070-independent/finite-prefix-area-receipt.json'
    fd, fr = json.loads(finite_result.read_text()), json.loads(finite_receipt.read_text())
    assert fr['exit_status'] == 0 and fr['owned_pid_absent'] and not fr['stderr']
    assert fr['result_sha256'] == sha(finite_result)
    assert fd['source_sha256'] == fr['source_sha256'] == sha(args.review_root/'computation/review_r070_finite_prefix_area.py')
    assert fd['input_sha256'] == fr['input_sha256'] == sha(args.review_root/'computation/inputs/r070-finite-prefix-area.json')
    result = dict(status='PASS_BYTE_AND_PROVENANCE_BINDINGS_ONLY', checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_sha256=sha(Path(__file__)), spec_sha256=sha(args.spec), packets=packet_results, frozen_report=spec['report'], recovered_runs=recovered,
                  new_finite_scalar_result_sha256=sha(finite_result),
                  new_finite_scalar_receipt_sha256=sha(finite_receipt),
                  saved_root_receipt_checks='Root prefix receipt: exit0,parent+8 result-returning worker PIDs recorded absent,12requested; closed pool described in source README. Wolfram receipt: exit0,PID absent,14saved checks all true. Historical receipt assertions read, not current PID observations.',
                  scope='Byte checks and provenance recovery only. Mathematical reading and independent computations are separately described in the review. No old run or PID/exit state was rerun/re-observed; no new halo or C2 claim.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], packets=len(packet_results), checked_packet_rows=sum(len(p['checked_files']) for p in packet_results), recovered_runs=len(recovered))))


if __name__ == '__main__':
    main()
