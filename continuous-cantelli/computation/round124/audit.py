from datetime import datetime
from fractions import Fraction as Q
from pathlib import Path
from statistics import median
import json
import re
import subprocess

import certificate as cert

ROOT=Path(__file__).parent


def mean_tent_upper():
    cert.exp_table()
    nx=512
    numerator=0
    for ix in range(-cert.R*nx,cert.R*nx):
        min_dist=max(0,ix-nx//2,nx//2-(ix+1))
        f_upper=max(Q(0),Q(7,4)-Q(min_dist,nx))
        f_num=int(f_upper*nx)
        min_x2,_=cert.square_bounds(ix,ix+1,0)
        numerator+=f_num*cert.exp_upper(min_x2,2*nx*nx)
    normal=cert.gamma_normalizer()[0]
    return Q(numerator,nx*nx*normal)


def review(report_file,metrics_file,kind):
    report=json.loads((ROOT/report_file).read_text())
    metrics=[json.loads(line) for line in (ROOT/metrics_file).read_text().splitlines()]
    pids=[metrics[0]['parent']]+metrics[0]['pids']
    alive=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    assert not alive.stdout.strip(),alive.stdout
    assert set(report['exits'].values())=={0}
    assert set(map(str,pids[1:]))==set(report['exits'])
    workers=report['stats'] if kind=='search' else report['workers']
    start_label,end_label=('started','ended') if kind=='search' else ('start_utc','end_utc')
    start=max(datetime.fromisoformat(r[start_label]) for r in workers)
    end=min(datetime.fromisoformat(r[end_label]) for r in workers)
    allr=[m for m in metrics[1:] if len(m['workers'])==12 and all(w.split()[2][0]=='R' for w in m['workers'])]
    cpu=[sum(float(w.split()[1]) for w in m['workers']) for m in allr]
    idle=[float(re.search(r'([\d.]+)% idle',m['load'][0])[1]) for m in allr]
    out=dict(parent=pids[0],worker_pids=pids[1:],worker_exit_codes=report['exits'],
             all_worker_window=[start.isoformat(),end.isoformat()],timed_seconds=(end-start).total_seconds(),
             all12_R_samples=len(allr),group_cpu_percent_median=median(cpu) if cpu else None,
             all12_R_system_idle_bounds=[min(idle),max(idle)] if idle else None,processes_after=[])
    if kind=='cert':
        data=report['raw']
        q=Q(data['integral_num'],report['nx']*data['weight_den']*data['normal_hi'])
        target=Q(data['target_hi'],data['scale']*data['weight_den'])
        tail=Q(data['tail'],data['scale']*data['weight_den'])
        delta=q-target-tail
        total_weight=sum(abs(Q(w)) for *_,w in report['tests'])
        assert delta==Q(report['certified_margin'])
        assert delta/total_weight==Q(report['max_residual_lower'])
        assert sum(r['integral_numerator'] for r in report['workers'])==data['integral_num']
        assert tail<Q(135,10**14)
        out|=dict(boxes=report['boxes'],exact_margin=str(delta),joint_residual_lower=str(delta/total_weight))
    return out


def main():
    mean=mean_tent_upper()
    assert mean<Q(23,25)
    scan=review('search-summary.json','telemetry.jsonl','search')
    cert1=review('certificate.json','certificate-telemetry.jsonl','cert')
    cert2=review('certificate-right.json','certificate-right-telemetry.jsonl','cert')
    assert Q(cert2['joint_residual_lower'])>Q(1,500)
    out=dict(search=scan,first_certificate=cert1,stronger_certificate=cert2,
             exact_tent_mean_upper=str(mean),tent_mean_upper_float=float(mean),ramp_mean_exact='5/4')
    (ROOT/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(tent_mean_upper=float(mean),sweep_seconds=scan['timed_seconds'],
                          certificate_seconds=[cert1['timed_seconds'],cert2['timed_seconds']])))


if __name__=='__main__':
    main()
