"""Local n16 Metal repeat comparison, with explicit wait and final record."""
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parent
old = root
binary = old / 'tau2_merge_rb_calibration_r026'
source = binary.with_suffix('.swift')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert sha(source) == '3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6'
raw = root / 'n16.jsonl'
assert not raw.exists()
domain = '9e52128000000001'
command = ['/usr/bin/time', '-l', str(binary), '16', '131072', '12', '2', '2',
           str(raw), domain, 'strat', '128']
state = dict(parameters=dict(n=16, P=131072, R=12, tau=2, L=128, domain=domain),
             source_sha256=sha(source), binary_sha256=sha(binary), command=command,
             started_utc=datetime.now(timezone.utc).isoformat())
with (root/'stdout.txt').open('x') as stdout, (root/'resource.txt').open('x') as resource:
    process = subprocess.Popen(command, stdout=stdout, stderr=resource)
    wake = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(process.pid)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    state.update(pid=process.pid, caffeinate_pid=wake.pid)
    (root/'started.json').write_text(json.dumps(state, indent=2)+'\n')
    status = process.wait()
    wake_status = wake.wait(timeout=60)
state.update(returncode=status, caffeinate_returncode=wake_status,
             waited_reaped=True, resource_sha256=sha(root/'resource.txt'))
if status == 0:
    rows = [json.loads(s) for s in raw.read_text().splitlines()]
    for kind in ('replicate', 'firstmerge_rb_calibration'):
        assert sorted(r['replicate'] for r in rows if r['kind']==kind)==list(range(12))
    assert sum(r['kind']=='run_end' for r in rows)==1
    assert rows[0]['domain']==domain and rows[0]['n']==16
    state.update(raw_sha256=sha(raw), rows=len(rows))
state['ended_utc']=datetime.now(timezone.utc).isoformat()
(root/'finished.json').write_text(json.dumps(state,indent=2)+'\n')
raise SystemExit(status)
