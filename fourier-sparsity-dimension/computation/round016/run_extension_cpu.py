#!/usr/bin/env python3
"""Run the R16 CPU verifier foreground, sampling only its owned child."""
import datetime as dt
import hashlib
import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def command(argv):
    p = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       check=False)
    return {'argv': argv, 'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}


def sample(pid, label):
    ps = command(['/bin/ps', '-p', str(pid), '-o', 'pid=,pcpu=,pmem=,etime=,command='])
    top = command(['/usr/bin/top', '-l', '1', '-n', '0'])
    cpu = next((line.strip() for line in top['stdout'].splitlines() if 'CPU usage:' in line), '')
    return {'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'label': label,
            'owned_child_ps': ps['stdout'].strip(), 'host_cpu_line': cpu}


def main():
    binary = ROOT / 'verify_extension_cpu'
    seed = ROOT / 'seed_masks.txt'
    output = ROOT / 'verify-result.json'
    receipt = ROOT / 'run-receipt.json'
    argv = [str(binary), str(seed), '15', str(output)]
    began = time.monotonic()
    record = {'schema': 'a28-r016-extension-cpu-run-v1',
              'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'command': argv,
              'hashes': {'runner_sha256': sha(__file__), 'cpp_sha256': sha(ROOT / 'verify_extension_cpu.cpp'),
                         'binary_sha256': sha(binary), 'seed_json_sha256': sha(ROOT / 'extension-seed.json'),
                         'seed_text_sha256': sha(seed), 'seed_extraction_sha256': sha(ROOT / 'seed-extraction.json')},
              'compiler': command(['/usr/bin/clang++', '--version']),
              'python': command([sys.executable, '--version']), 'samples': []}
    proc = subprocess.Popen(argv, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    record['child_pid'] = proc.pid
    record['samples'].append(sample(proc.pid, 'start'))
    next_sample = time.monotonic() + 0.5
    while proc.poll() is None:
        if time.monotonic() - began > 600:
            proc.terminate()
            out, err = proc.communicate(timeout=15)
            record['timed_out'] = True
            break
        if time.monotonic() >= next_sample:
            record['samples'].append(sample(proc.pid, 'interval'))
            next_sample += 0.5
        time.sleep(0.05)
    else:
        out, err = proc.communicate()
        record['timed_out'] = False
    record['samples'].append(sample(proc.pid, 'end'))
    record.update({'finished_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                   'elapsed_seconds_wall': time.monotonic() - began,
                   'child_returncode': proc.returncode, 'child_stdout': out,
                   'child_stderr': err,
                   'result_sha256': sha(output) if output.exists() else None})
    if output.exists():
        record['result'] = json.loads(output.read_text())
    receipt.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    return 0 if proc.returncode == 0 and not record['timed_out'] else (proc.returncode or 20)


if __name__ == '__main__':
    sys.exit(main())
