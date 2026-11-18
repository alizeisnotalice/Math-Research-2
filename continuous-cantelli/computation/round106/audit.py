from datetime import datetime
from math import log, pi
from pathlib import Path
from statistics import median
import json
import re
import subprocess

from scipy.optimize import brentq

from multiscale_mellin import evaluate

ROOT=Path(__file__).parent


def lock_phase(l0,beta,omega,theta):
    central=l0-log(2*l0)
    turns=round((omega*central**beta-theta)/(2*pi))
    ell=((2*pi*turns+theta)/omega)**(1/beta)
    return brentq(lambda x:x-log(2*x)-ell,ell,ell+2*log(2*ell))


def main():
    s=json.loads((ROOT/'summary.json').read_text())
    rows=[]
    for beta in [.8,1.,1.25,1.5,1.8]:
        for l0 in [64,1024,8192,65536,262144]:
            for theta in [0.,pi]:
                if beta>1.5 and l0>65536:
                    continue
                omega,amp=.8,.8
                l=lock_phase(l0,beta,omega,theta)
                row=evaluate(-1,(l,beta,omega,theta,amp),check=True)
                rows.append(row)
    samples=[json.loads(x) for x in (ROOT/'telemetry.jsonl').read_text().splitlines()[1:]]
    all12=[line for line in samples if len(line['ps'])==12 and all(r.split()[2][0]=='R' for r in line['ps'])]
    group=[sum(float(r.split()[1]) for r in line['ps']) for line in all12]
    idle=[float(re.search(r'([\d.]+)% idle',line['load'][0])[1]) for line in all12]
    latest=max(datetime.fromisoformat(w['start']) for w in s['workers'])
    earliest=min(datetime.fromisoformat(w['end']) for w in s['workers'])
    pids=[s['parent']]+s['pids']
    query=subprocess.run(['ps','-p',','.join(str(pid) for pid in pids),'-o','pid='],
                         capture_output=True,text=True)
    out=dict(phase_locked=rows,phase_parameters='one fixed core for each beta, omega=amp=.8, phi=0',
             all12_running_samples=len(all12),idle_range=[min(idle),max(idle)],
             cpu_percent_median=median(group),cpu_percent_range=[min(group),max(group)],
             all12_work_window=[latest.isoformat(),earliest.isoformat()],
             all12_work_seconds=(earliest-latest).total_seconds(),processes_after=query.stdout.splitlines())
    (ROOT/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='phase_locked'}))


if __name__=='__main__':
    main()
