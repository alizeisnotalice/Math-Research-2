#!/usr/bin/env python3
"""Foreground ownership/receipt wrapper for the bounded R19 CPU scan."""
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(argv):
    p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       text=True, check=False)
    return {'argv': argv, 'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}


def sample(pid, label):
    ps = run(['/bin/ps', '-p', str(pid), '-o', 'pid=,pcpu=,pmem=,etime=,command='])
    top = run(['/usr/bin/top', '-l', '1', '-n', '0'])
    line = next((z.strip() for z in top['stdout'].splitlines() if 'CPU usage:' in z), '')
    mem = next((z.strip() for z in top['stdout'].splitlines() if 'PhysMem:' in z), '')
    return {'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'label': label,
            'owned_child_ps': ps['stdout'].strip(), 'host_cpu_line': line, 'host_memory_line': mem}


def main():
    binary = ROOT / 'audit_rm36_weight16'
    output = ROOT / 'full-result.json'
    receipt = ROOT / 'full-run-receipt.json'
    argv = [str(binary), 'full', '15', str(output)]
    record = {'schema': 'a28-r019-rm36-full-run-v1', 'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'command': argv, 'threads': 15,
              'hashes': {'source_sha256': sha(ROOT/'audit_rm36_weight16.cpp'), 'binary_sha256': sha(binary),
                         'runner_sha256': sha(__file__), 'orbit_input_sha256': sha(ROOT/'orbits-gpu.json')},
              'compiler': run(['/usr/bin/clang++', '--version']), 'python': run([sys.executable, '--version']),
              'samples': []}
    start = time.monotonic()
    p = subprocess.Popen(argv, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    record['child_pid'] = p.pid
    record['samples'].append(sample(p.pid, 'start'))
    next_sample = time.monotonic() + 1
    while p.poll() is None:
        if time.monotonic() - start > 900:
            p.terminate(); out, err = p.communicate(timeout=15); record['timed_out'] = True; break
        if time.monotonic() >= next_sample:
            record['samples'].append(sample(p.pid, 'interval')); next_sample += 1
        time.sleep(0.05)
    else:
        out, err = p.communicate(); record['timed_out'] = False
    record['samples'].append(sample(p.pid, 'end'))
    record.update({'finished_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                   'elapsed_seconds_wall': time.monotonic()-start, 'child_returncode': p.returncode,
                   'child_stdout': out, 'child_stderr': err,
                   'result_sha256': sha(output) if output.exists() else None})
    if output.exists(): record['result'] = json.loads(output.read_text())
    receipt.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    return 0 if p.returncode == 0 and not record['timed_out'] else (p.returncode or 20)


if __name__ == '__main__':
    sys.exit(main())
