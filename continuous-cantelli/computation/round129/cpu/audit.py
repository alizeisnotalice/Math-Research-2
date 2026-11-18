from datetime import datetime
from fractions import Fraction as Q
from pathlib import Path
import json
import subprocess

import bound_primitives as b
import refine

ROOT=Path(__file__).parent


def verify(name):
    a=json.loads((ROOT/f'{name}.json').read_text())
    b.EXP_GRID=a.get('exponential_grid',4096)
    b.EXP_END=40*b.EXP_GRID
    b.exp_table()
    assert len(a['workers'])==12
    assert set(a['exits'].values())=={0}
    assert set(a['exits'])==set(map(str,a['pids']))
    pids=[a['parent']]+a['pids']
    ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid='],capture_output=True,text=True)
    assert not ps.stdout.splitlines()
    lo,hi=b.gamma_normalizer()
    tl,tu,tn=b.targets_and_tail()
    tp=refine.positive_tail()
    den=b.SCALE*b.WEIGHT_DEN
    l1=Q(sum(abs(w) for *_,w in b.TESTS),b.WEIGHT_DEN)
    for i,r in enumerate(a['brackets']):
        lo_raw=sum(w['lower'][i] for w in a['workers'])
        hi_raw=sum(w['upper'][i] for w in a['workers'])
        left=Q(lo_raw,a['nx']*b.WEIGHT_DEN*hi)-Q(tu+tn,den)
        right=Q(hi_raw,a['nx']*b.WEIGHT_DEN*lo)-Q(tl,den)+tp
        assert left==Q(r['lower'])<=Q(r['upper'])==right
        r['residual_lower']=str(left/l1)
        r['residual_lower_display']=float(left/l1)
        r['amplitude_if_left_ge_3_4']=str(Q(r['a'])-Q(3,4))
    last_start=max(datetime.fromisoformat(w['started']) for w in a['workers'])
    first_finish=min(datetime.fromisoformat(w['ended']) for w in a['workers'])
    return dict(worker_pids=a['pids'],parent=a['parent'],worker_exits=a['exits'],processes_after=[],
                gates=[last_start.isoformat(),first_finish.isoformat()],
                full_pool_elapsed=(first_finish-last_start).total_seconds(),
                group_cpu_median=a['group_cpu_median'],system_idle=a['host_idle'],stats=a['stats'],
                bounds=a['brackets'],tail_negative=float(Q(tn,den)),tail_positive=float(tp))


def main():
    out={k:verify(k) for k in ('target','fine')}
    exact={row['a']:row for row in out['fine']['bounds']}
    assert Q(exact['171/100']['lower'])>Q(1256,10**6)
    assert Q(exact['3417/2000']['lower'])>Q(429,10**6)
    (ROOT/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:{'window':v['gates'],'seconds':v['full_pool_elapsed'],'boxes':v['stats']['boxes']}
                      for k,v in out.items()}))


if __name__=='__main__':
    main()
