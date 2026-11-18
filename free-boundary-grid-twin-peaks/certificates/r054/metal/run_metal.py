"""Run the transport diagnostic with per-process and per-file receipts."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
from argparse import ArgumentParser
from subprocess import Popen,DEVNULL
import json

ROOT=Path(__file__).resolve().parent
SETTINGS={
    'gate': (3,17,1,0,'9e53030000000001'),
    'replay':(16,131072,12,16,'9e52128000000001'),
    'contact-gate':(3,17,1,0,'9e53030000000001'),
    'contact-n16':(16,131072,12,16,'9e52128000000001'),
    'contact-n24':(24,65536,8,16,'9e54100000000001'),
    'n24':(24,65536,8,64,'9e53640000000001')
}

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def main():
    p=ArgumentParser()
    p.add_argument('mode',choices=SETTINGS)
    args=p.parse_args()
    n,P,R,L,domain=SETTINGS[args.mode]
    out=ROOT/args.mode
    out.mkdir(exist_ok=True)
    raw=out/'run.jsonl'
    assert not raw.exists()
    name='contact' if args.mode.startswith('contact') else 'transport'
    binary=ROOT/name
    source=ROOT/(name+'.swift')
    command=['/usr/bin/time','-l',str(binary),str(n),str(P),str(R),'2','2',
             str(raw),domain,'strat',str(L)]
    if args.mode in ('gate','contact-gate'):
        command+=['trace']
    state={'mode':args.mode,'command':command,'parameters':dict(n=n,P=P,R=R,tau=2,L=L,domain=domain),
           'source_sha256':digest(source),'binary_sha256':digest(binary),'started_utc':now()}
    with (out/'stdout.txt').open('x') as stdout, (out/'resource.txt').open('x') as resource:
        with Popen(command,stdout=stdout,stderr=resource,cwd=out) as proc:
            wake=Popen(['/usr/bin/caffeinate','-i','-w',str(proc.pid)],stdout=DEVNULL,stderr=DEVNULL)
            state.update(pid=proc.pid,caffeinate_pid=wake.pid)
            (out/'started.json').write_text(json.dumps(state,indent=2)+'\n')
            status=proc.wait()
            wake_status=wake.wait(timeout=60)
    state.update(ended_utc=now(),returncode=status,caffeinate_returncode=wake_status,
                 waited_reaped=True,resource_sha256=digest(out/'resource.txt'))
    if status==0:
        payload=[json.loads(x) for x in raw.read_text().splitlines()]
        assert payload[-1]['kind']=='run_end' and len(payload)==2*R+2
        assert all(len(x['transport_histograms'])==4 for x in payload[2:-1:2])
        state.update(raw_sha256=digest(raw),rows=len(payload))
    (out/'finished.json').write_text(json.dumps(state,indent=2)+'\n')
    raise SystemExit(status or wake_status)

if __name__=='__main__':
    main()
