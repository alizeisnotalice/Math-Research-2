#!/usr/bin/env python3
"""Replay all 11,760 corrected S={2} cores and both trees from empty storage."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
from common import ROOT, canonical, guard_assertions, require, sha
from verify_certificate import audit_run, verify_manifest


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def main():
    guard_assertions()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path, help='new directory outside the certificate')
    parser.add_argument('--workers', type=int, choices=range(1, 14), default=4)
    options = parser.parse_args()
    out = options.out.resolve()
    require(out != ROOT and ROOT not in out.parents, '--out must be outside the sealed certificate')
    if options.out.exists() or options.out.is_symlink():
        raise SystemExit('refuse existing --out; a fresh run requires a new directory')
    import sympy
    manifest_hash = verify_manifest(ROOT)
    out.mkdir(parents=True, exist_ok=False)
    scratch = out / 'scratch'
    scratch.mkdir()
    logs = out / 'logs'
    logs.mkdir()
    stage_runs = []
    helpers = ROOT / 'modules'
    environment = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'TMPDIR': str(scratch),
                   'TMP': str(scratch), 'TEMP': str(scratch)}
    began, utc_began = time.monotonic(), utc_now()

    def run(label, module, *args):
        program = helpers / module
        command = [sys.executable, str(program), *map(str, args)]
        started = time.monotonic()
        stdout, stderr = logs / f'{label}.stdout.txt', logs / f'{label}.stderr.txt'
        print(f'{utc_now()} {label}', flush=True)
        with stdout.open('xb') as output_stream, stderr.open('xb') as error_stream:
            with subprocess.Popen(command, env=environment, cwd=out, stdout=output_stream, stderr=error_stream) as process:
                pid = process.pid
                exit_code = process.wait()
        item = {'stage': label, 'command': command, 'code_sha256': sha(program), 'pid': pid, 'exit_code': exit_code,
                'seconds': time.monotonic()-started, 'stdout_sha256': sha(stdout), 'stderr_sha256': sha(stderr)}
        (logs / f'{label}.json').write_bytes(canonical(item))
        stage_runs.append(item)
        if exit_code:
            raise subprocess.CalledProcessError(exit_code, command)

    run('stage4', 's2_degree8_stage4.py',
        '--out', out / 'stage4', '--workers', options.workers, '--max-seconds', 14400)
    run('stage58', 's2_degree8_stage5to8.py', '--stage4', out / 'stage4',
        '--out', out / 'stage58', '--workers', options.workers, '--max-seconds', 14400)
    run('root', 's2_degree8_root_screen.py', '--stage58', out / 'stage58',
        '--out', out / 'root', '--workers', options.workers, '--max-seconds', 14400)
    run('tree-input', 'extract_tree_input.py', '--root', out / 'root', '--out', out / 'tree-input.jsonl.gz')
    run('trees', 's2_degree8_dual_trees.py', '--input', out / 'tree-input.jsonl.gz',
        '--out', out / 'trees', '--workers', options.workers, '--depth', 128, '--cap', 25000, '--tree-seconds', 20)

    original = audit_run()
    replica = audit_run(out)
    require(original['semantics'] == replica['semantics'], 'fresh replay differs from complete sealed records')
    require(manifest_hash == verify_manifest(ROOT), 'certificate changed during fresh replay')
    report = {'status': 'passed_fresh_full_replay', 'package_manifest_sha256': manifest_hash,
              'cached_cores': 0, 'workers': options.workers, 'seed': None, 'numeric_tolerance': None,
              'python': sys.version, 'python_executable': sys.executable, 'sympy': sympy.__version__,
              'started_at': utc_began, 'finished_at': utc_now(), 'seconds': time.monotonic()-began,
              'commands': stage_runs, 'original': original, 'replay': replica}
    scratch.rmdir()
    (out / 'fresh-report.json').write_bytes(canonical(report))
    print(json.dumps({'status': report['status'], 'seconds': report['seconds'], 'trees': replica['semantics']['trees']}, sort_keys=True))


if __name__ == '__main__':
    main()
