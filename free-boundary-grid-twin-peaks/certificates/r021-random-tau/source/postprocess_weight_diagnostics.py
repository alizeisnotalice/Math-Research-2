#!/usr/bin/env python3
"""Manifest-replayed log-weight concentration diagnostics for archived r021 runs.

This deliberately does not alter raw output or make a probability claim.
"""
import argparse, json, math, pathlib, sys
BASE=0xA022000000000000; MASK=(1<<64)-1
def lse(xs):
    m=max(xs); return m+math.log(sum(math.exp(x-m) for x in xs))
def diagnostic(logw):
    a=lse(logw); b=lse([2*x for x in logw])
    return {'weighted_ESS':math.exp(2*a-b),'max_normalized_weight':math.exp(max(logw)-a),
            'min_normalized_weight':math.exp(min(logw)-a)}
def main():
    p=argparse.ArgumentParser();p.add_argument('raw');p.add_argument('--n',type=int,required=True);p.add_argument('--K',type=int,default=6);a=p.parse_args()
    raw=pathlib.Path(a.raw);manifest=json.load(open(str(raw)+'.manifest.json'))
    if manifest.get('n')!=a.n or not all(k in manifest for k in ('reps','dispatch')):raise RuntimeError('manifest mismatch')
    here=str(pathlib.Path(__file__).parent)
    if here not in sys.path:sys.path.insert(0,here)
    import random_tau_design as d
    ds=manifest['dispatch']; R=manifest['reps']
    if len(ds)!=R or sorted(x.get('globalrep',-1) for x in ds)!=list(range(R)):raise RuntimeError('manifest globalrep')
    td={x.get('tau_domain_hex') for x in ds};sd={x.get('trajectory_domain_hex') for x in ds}
    if len(td)!=1 or len(sd)!=1:raise RuntimeError('manifest domains')
    expected=d.dispatch_records(a.n,R,a.K,int(next(iter(td)),16),int(next(iter(sd)),16))[1]
    ek=('tau','qnum','qden','tau_u','tau_rejections','desired_seed_hex','domain_hex','tau_domain_hex','trajectory_domain_hex')
    eb={x['globalrep']:x for x in expected};mb={x['globalrep']:x for x in ds}
    if any(any(mb[i].get(k)!=eb[i][k] for k in ek) for i in eb):raise RuntimeError('manifest replay')
    rows=[json.loads(x) for x in open(raw) if x.strip()]
    if len(rows)!=R or {x.get('globalrep') for x in rows}!=set(eb):raise RuntimeError('raw incomplete')
    if len({x.get('desired_seed_hex') for x in rows})!=R:raise RuntimeError('duplicate desired seed')
    N=a.n*a.n; wz=[];wd=[];dead=0
    for x in rows:
        if any(x.get(k)!=eb[x['globalrep']][k] for k in ek):raise RuntimeError('row replay')
        if x.get('seed_hex','').lower().zfill(16)!=x['desired_seed_hex'].lower().zfill(16):raise RuntimeError('seed mismatch')
        if (BASE+int(x['domain_hex'],16)+a.n*1_000_000+x['tau']*10_000)&MASK!=int(x['desired_seed_hex'],16):raise RuntimeError('CLI inverse')
        if x['dead']:
            dead+=1;continue
        q=x['qnum']/x['qden'];lz=float(x['logZ'])-math.log(q)
        wz.append(lz)
        u=float(x['U_over_Z'])
        if u>0:wd.append(lz+math.log(N*u))
    if not wz:raise RuntimeError('all observations dead')
    out={'kind':'importance_weight_diagnostic','n':a.n,'observations':R,'nondead':len(wz),'dead':dead,
         'Z_weight':diagnostic(wz),'D_weight':None if not wd else diagnostic(wd),
         'note':'ESS/max weight diagnose finite importance-weight concentration only. They are not independent-particle counts, error bars, or conditional-law probabilities.'}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
