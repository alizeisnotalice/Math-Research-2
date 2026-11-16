"""Foreground child runner; telemetry is device-wide, not exclusive GPU usage."""
import json
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

output = Path(sys.argv[1])
command = sys.argv[2:]
assert not output.exists(), 'Never overwrite an existing run telemetry file'
process = subprocess.Popen(command)
stop = threading.Event()

def sample():
    data = {'kind': 'sample', 'epoch': time.time(), 'pid': process.pid}
    try:
        data['rss_kib'] = int(subprocess.check_output(
            ['ps', '-o', 'rss=', '-p', str(process.pid)], text=True).strip())
    except (subprocess.CalledProcessError, ValueError):
        data['rss_kib'] = None
    try:
        raw = subprocess.check_output(['ioreg', '-r', '-c', 'AGXAccelerator', '-d', '1'], text=True)
        for field in ['Device Utilization %', 'Renderer Utilization %', 'Alloc system memory', 'In use system memory']:
            match = re.search('"' + re.escape(field) + r'"=(\d+)', raw)
            data[field] = int(match.group(1)) if match else None
    except subprocess.CalledProcessError:
        data['gpu_telemetry_error'] = True
    return data

with output.open('x') as handle:
    handle.write(json.dumps({'kind': 'start', 'pid': process.pid, 'command': command,
                             'period_seconds': 5, 'scope': 'target RSS; whole-device AGX, nonexclusive'}) + '\n')
    handle.flush()
    def loop():
        while not stop.is_set():
            handle.write(json.dumps(sample()) + '\n')
            handle.flush()
            stop.wait(5)
    thread = threading.Thread(target=loop)
    thread.start()
    returncode = process.wait()
    stop.set()
    thread.join()
    handle.write(json.dumps({'kind': 'exit', 'pid': process.pid, 'returncode': returncode,
                             'epoch': time.time()}) + '\n')
raise SystemExit(returncode)
