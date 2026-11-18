from collections import Counter
from datetime import datetime
from pathlib import Path
from statistics import median
import json
import re
import subprocess

import mpmath as mp

from bounded_heat import Image, cap_opt, ndtr, phi

ROOT=Path(__file__).parent


def high_precision(row):
    mp.mp.dps=55
    a,b,t=map(lambda k:mp.mpf(str(row[k])),('a','b','t'))
    mass=(mp.erfc(a/mp.sqrt(2))+mp.erfc(b/mp.sqrt(2)))/2
    normal=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
    mean=normal(b)-normal(a)
    width=a+b
    wave=[n*mp.pi/width for n in range(1,1+int(mp.ceil(mp.sqrt(270*width**2/(mp.pi**2*t)))))]
    coeff=[2/width*mp.sin(w*a)*mp.exp(-t*w*w/2) for w in wave]

    def primitive(x,j):
        values=[]
        for w,c_i in zip(wave,coeff):
            c,s=mp.cos(w*(x+a)),mp.sin(w*(x+a))
            if j==0:
                term=-c/w
            elif j==1:
                term=-x*c/w+s/(w*w)
            else:
                term=-x*x*c/w+2*x*s/(w*w)+2*c/(w**3)
            values.append(c_i*term)
        return mp.fsum(values)

    base_l=[primitive(-a,j) for j in range(3)]
    base_r=[primitive(b,j) for j in range(3)]

    def two_tail(l,u,j):
        return primitive(l,j)-base_l[j]+base_r[j]-primitive(u,j)

    l,u=mp.findroot(lambda l,u:(two_tail(l,u,0)-mass,two_tail(l,u,1)-mean),
                    tuple(mp.mpf(str(x)) for x in row['cuts']))
    exact_gain=a*b*mass+(b-a)*mean-two_tail(l,u,2)
    reported=mp.mpf(str(row['gain']))
    assert abs(exact_gain-reported)<mp.mpf('3e-10'), (row,exact_gain)
    return dict(a=str(a),b=str(b),t=str(t),
                gain=mp.nstr(exact_gain,46),cuts=[mp.nstr(l,32),mp.nstr(u,32)],
                discrepancy=mp.nstr(exact_gain-reported,12))


def main():
    rows=[row for path in ROOT.glob('shard-*.json') for row in json.loads(path.read_text())['rows']]
    feasible=[row for row in rows if row['reason']=='feasible']
    counts=Counter()
    for row in rows:
        support=(-row['a']*row['mass']<=row['mean']<=row['b']*row['mass'])
        counts[(row['reason'],support)]+=1
    def available(row):
        return 1-(1+row['t'])*row['mass']

    previous=[row for row in feasible if row['base']>available(row)]
    new=[row for row in feasible if row['base']<=available(row)<row['base']+row['gain']]
    biggest=max(new,key=lambda row:row['base']+row['gain']-available(row),default=None)
    feasible.sort(key=lambda row:row['relative_gain'])
    percentiles=[0,.1,.5,.9,1]
    targets=[feasible[round(x*(len(feasible)-1))] for x in percentiles]
    demos=[]
    for a,b,t in [(.5,.5,.05),(1.,1.,.25),(1.,1.,.75),(2.,2.,.5),
                  (.5,1.5,.15),(.6,1.5,.15),(.8,1.5,.25)]:
        mass=ndtr(-a)+ndtr(-b)
        mean=phi(b)-phi(a)
        base=a*(phi(b)-b*ndtr(-b))+b*(phi(a)-a*ndtr(-a))
        val=cap_opt(Image(a,b,t),mass,mean)|dict(a=a,b=b,t=t,mass=float(mass),mean=float(mean),base=float(base))
        if val['reason']=='feasible':
            val['relative_gain']=val['gain']/base
        demos.append(val)
    summary=json.loads((ROOT/'summary.json').read_text())
    metrics=[json.loads(line) for line in (ROOT/'telemetry.jsonl').read_text().splitlines()][1:]
    all12=[m for m in metrics if len(m['workers'])==12 and all(r.split()[2][0]=='R' for r in m['workers'])]
    usage=[sum(float(w.split()[1]) for w in m['workers']) for m in all12]
    idle=[float(re.search(r'([\d.]+)% idle',m['load'][0])[1]) for m in all12]
    start=max(datetime.fromisoformat(s['start_utc']) for s in summary['shards'])
    end=min(datetime.fromisoformat(s['end_utc']) for s in summary['shards'])
    all_pids=[summary['parent']]+summary['pids']
    ps=subprocess.run(['ps','-p',','.join(str(pid) for pid in all_pids),'-o','pid='],
                      capture_output=True,text=True)
    result=dict(coverage={f'{k[0]}-support-{"ok" if k[1] else "fails"}':v for k,v in sorted(counts.items())},
                width_range=[min(r['a']+r['b'] for r in rows),max(r['a']+r['b'] for r in rows)],
                time_range=[min(r['t'] for r in rows),max(r['t'] for r in rows)],
                samples=targets,demos=demos,high_precision=[high_precision(row) for row in targets],
                normalization_check=dict(base_alone_excluded=len(previous),
                                         additional_with_cap=len(new),largest_new_example=biggest),
                telemetry=dict(all_12_R_samples=len(all12),idle_range=[min(idle),max(idle)],
                               group_cpu_median=median(usage),group_cpu_max=max(usage),
                               all_worker_timed_interval=[start.isoformat(),end.isoformat()],
                               interval_seconds=(end-start).total_seconds(),processes_after=ps.stdout.splitlines()))
    (ROOT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(coverage=result['coverage'],demos=result['demos'],telemetry=result['telemetry'],
                          high_precision=result['high_precision'])))


if __name__=='__main__':
    main()
