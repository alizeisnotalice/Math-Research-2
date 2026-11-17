#!/usr/bin/env python3
"""Replay finite audits and small C++ gates; never rerun n7 production."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--module-cache', type=Path, required=True)
    p.add_argument('--compiler', default='c++')
    a = p.parse_args()
    root = Path(__file__).resolve().parent
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    a.module_cache.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
               CLANG_MODULE_CACHE_PATH=str(a.module_cache.resolve()))
    records = []
    def run(label, cmd, cwd=out):
        began = time.time()
        with (out/(label+'.stdout')).open('w') as f, (out/(label+'.stderr')).open('w') as e:
            proc = subprocess.Popen(list(map(str, cmd)), cwd=cwd, env=env, stdout=f, stderr=e)
            rc = proc.wait()
        absent = subprocess.run(['ps', '-p', str(proc.pid), '-o', 'pid='],
                                capture_output=True).returncode != 0
        records.append(dict(label=label, pid=proc.pid, returncode=rc,
                            waited_reaped=True, absent_after=absent, seconds=time.time()-began))
        (out/'receipt.json').write_text(json.dumps(dict(python=sys.version,
                    records=records), indent=2)+'\n')
        assert rc == 0 and absent, (label, rc)
    for group, name in [('gap-mixture','audit_fixed_gap_mixture.py'),
                        ('weighted-reinsertion','audit_psi_T_n3.py'),
                        ('completion-transfer','audit_T_ge_S_n3.py'),
                        ('joint-geometry','audit_kernel.py')]:
        folder = out/group
        folder.mkdir()
        shutil.copy2(root/group/name, folder/name)
        run(group, [sys.executable, '-B', folder/name], folder)
    spatial = root/'spatial'
    run('n7-crt', [sys.executable, '-B', spatial/'scripts/audit_spatial_crt.py',
                   '--root', spatial, '--output', out/'n7-spatial-crt-audit.json'])
    assert json.loads((out/'n7-spatial-crt-audit.json').read_text()) == json.loads(
        (spatial/'evidence/n7-spatial-crt-audit.json').read_text())
    run('n7-residuals', [sys.executable, '-B', spatial/'scripts/audit_spatial_residuals.py',
                        '--root', spatial, '--output', out/'n7-spatial-residuals-audit.json'])
    assert json.loads((out/'n7-spatial-residuals-audit.json').read_text()) == json.loads(
        (spatial/'evidence/n7-spatial-residuals-audit.json').read_text())
    source = spatial/'source/gap_spatial_kernel_mod.cpp'
    run('compile', [a.compiler, '-O3', '-std=c++20', '-pthread',
                    '-fmodules-cache-path='+str(a.module_cache.resolve()), source,
                    '-o', out/'spatial-kernel'])
    gate = spatial/'spatial-small-gate'
    manifest = gate/'source/GATE_INPUTS.json'
    config = json.loads(manifest.read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == config['source_sha256']
    raw = out/'small-gate-raw'
    raw.mkdir()
    for n, vertices in config['gate_points'].items():
        for vertex in vertices:
            for pi in range(4):
                name = f'n{n}-v{vertex}-p{pi}'
                run(name, [out/'spatial-kernel', n, vertex, pi, 1, raw/(name+'.json')])
    run('small-gate-exact', [sys.executable, '-B', gate/'source/exact_spatial_gate.py',
                            '--runs', raw, '--references', gate/'references',
                            '--source', source, '--manifest', manifest,
                            '--output', out/'small-gate-exact.json'])
    print(json.dumps(dict(status='passed', processes=len(records),
                          n7_production_rerun=False, small_residues=60)))

if __name__ == '__main__':
    main()
