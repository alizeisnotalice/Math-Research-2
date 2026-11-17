#!/usr/bin/env python3
"""Detached launcher for the reviewed four-core, four-hour local54 complement run."""
from __future__ import annotations
import hashlib, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round047-local54-complement')
PRODUCER = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/local54_complement47/local54_complement.py')

def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

if not Path('/Volumes/data').is_mount() or not (ROOT/'input-manifest.json').is_file():
    raise SystemExit('data or prepared local54 manifest unavailable')
if any('local54_complement.py run' in x for x in subprocess.check_output(['ps','-axo','command='], text=True).splitlines()):
    raise SystemExit('local54 producer already live')
log = (ROOT/'run'/'local54-producer.log').open('a', buffering=1)
cmd = [sys.executable, str(PRODUCER), 'run', '--root', str(ROOT), '--workers', '4', '--seconds', '14400']
p = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                     start_new_session=True, close_fds=True)
time.sleep(2)
if p.poll() is not None:
    raise SystemExit(f'producer exited during launch: rc={p.returncode}')
out = {'schema':1, 'pid':p.pid, 'ppid_at_launch':os.getpid(), 'command':cmd, 'start_new_session':True,
       'producer_sha256':sha(PRODUCER), 'input_manifest_sha256':sha(ROOT/'input-manifest.json'), 'started_unix':time.time()}
(ROOT/'run'/'local54-launch.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
