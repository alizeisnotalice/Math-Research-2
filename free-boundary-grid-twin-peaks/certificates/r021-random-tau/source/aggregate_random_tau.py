#!/usr/bin/env python3
"""Strict Decimal-80 aggregator for random-tau r021 JSONL output."""
import argparse, json, pathlib, sys
from decimal import Decimal, getcontext
getcontext().prec=80
MASK=(1<<64)-1; GPU_BASE=0xA022000000000000
def D(x): return Decimal(str(x))
def fail(s): raise RuntimeError(s)
def finite(x):
    try: return D(x).is_finite()
    except Exception: return False
def recover(logz, scaled, factor=1):
    if logz is None or scaled is None or D(scaled)==0: return Decimal(0)
    return (D(logz)+D(scaled).ln()+D(factor).ln()).exp()
def svar(xs):
    if len(xs)<2:return None
    m=sum(xs)/len(xs);return sum((x-m)**2 for x in xs)/(len(xs)-1)
def main():
    p=argparse.ArgumentParser();p.add_argument('jsonl');p.add_argument('--n',type=int,required=True);p.add_argument('--K',type=int,default=6);p.add_argument('--manifest');a=p.parse_args()
    mp=pathlib.Path(a.manifest or (a.jsonl+'.manifest.json'))
    if not mp.exists():fail('missing manifest')
    m=json.load(open(mp)); N=a.n*a.n;L=N-1;K=min(a.K,N-2);tail=N-K-1;qden=4*(1<<K)*tail*L
    if m.get('n')!=a.n or not all(k in m for k in ('P','moveL','method','reps','dispatch')):fail('bad manifest run fields')
    ds=m['dispatch']; R=m['reps']
    if R<=0 or len(ds)!=R or sorted(d.get('globalrep',-1) for d in ds)!=list(range(R)):fail('bad manifest globalrep set')
    tau_domains={d.get('tau_domain_hex') for d in ds};seed_domains={d.get('trajectory_domain_hex') for d in ds}
    if len(tau_domains)!=1 or len(seed_domains)!=1:fail('manifest domains not constant')
    # Recreate its complete stochastic dispatch from independent namespaces.
    here=str(pathlib.Path(__file__).with_name('random_tau_design.py').parent)
    if here not in sys.path:sys.path.insert(0,here)
    import random_tau_design as design
    replay=design.dispatch_records(a.n,R,K,int(next(iter(tau_domains)),16),int(next(iter(seed_domains)),16))[1]
    rb={x['globalrep']:x for x in replay};gb={x['globalrep']:x for x in ds}
    dk=('tau','qnum','qden','tau_u','tau_rejections','desired_seed_hex','domain_hex','tau_domain_hex','trajectory_domain_hex')
    if any(any(gb[i].get(k)!=rb[i][k] for k in dk) for i in rb):fail('manifest dispatch does not replay')
    seen=set();seeds=set();rows=[]
    for line in open(a.jsonl):
        r=json.loads(line)
        if r.get('kind')!='replicate':continue
        need=('globalrep','tau','qnum','qden','tau_u','tau_rejections','desired_seed_hex','domain_hex','seed_hex','logZ','dead','n','P','moveL','method','Z','U_over_Z','scaled_U','D','F','G','F_over_Z','G_over_Z','F_minus_2G_minus_4D')
        if any(k not in r for k in need):fail('row missing required field')
        gr=r['globalrep']
        if gr in seen:fail('duplicate globalrep')
        if r['desired_seed_hex'] in seeds:fail('duplicate desired seed')
        if gr not in rb:fail('unexpected globalrep')
        seen.add(gr);seeds.add(r['desired_seed_hex'])
        if any(r.get(k)!=rb[gr][k] for k in dk):fail('row dispatch mismatch')
        if (r['n'],r['P'],r['moveL'],r['method'])!=(a.n,m['P'],m['moveL'],m['method']):fail('row run metadata mismatch')
        if r['seed_hex'].lower().zfill(16)!=r['desired_seed_hex'].lower().zfill(16):fail('seed mismatch')
        j=r['tau']-1;num=3*L*((1<<(K-j))*tail if j<=K else 1)+(1<<K)*tail
        if not (2<=r['tau']<=N and r['qnum']==num and r['qden']==qden):fail('wrong q metadata')
        if (GPU_BASE+int(r['domain_hex'],16)+a.n*1_000_000+r['tau']*10_000)&MASK != int(r['desired_seed_hex'],16):fail('CLI domain mismatch')
        q=D(r['qnum'])/D(r['qden'])
        if r['dead']:
            if r['Z']!=0 or any(r[k] is not None for k in ('logZ','U_over_Z','scaled_U','D','F','G','F_over_Z','G_over_Z','F_minus_2G_minus_4D')):fail('dead payload is not canonical zero/null')
            z=dist=f=g=Decimal(0)
        else:
            numeric=('logZ','U_over_Z','scaled_U','D','F','G','F_over_Z','G_over_Z','F_minus_2G_minus_4D')
            if any(r[k] is None or not finite(r[k]) for k in numeric):fail('nondead null/nonfinite field')
            if D(r['scaled_U'])!=D(r['U_over_Z']):fail('scaled_U mismatch')
            z=D(r['logZ']).exp()
            if r['Z']!=0 and abs(D(r['Z'])-z)>D('1e-12')*max(abs(D(r['Z'])),abs(z)):fail('Z/logZ mismatch')
            dist=recover(r['logZ'],r['U_over_Z'],N);f=recover(r['logZ'],r['F_over_Z']);g=z*D(r['G_over_Z'])
            resid=D(r['F_over_Z'])-2*D(r['G_over_Z'])-4*D(N)*D(r['U_over_Z'])
            scale=max(D(1),abs(D(r['F_over_Z'])),abs(2*D(r['G_over_Z'])),abs(4*D(N)*D(r['U_over_Z'])))
            if abs(resid)>D('1e-8')*scale:fail('F-2G-4D scaled identity')
            rawres=D(r['F'])-2*D(r['G'])-4*D(r['D'])
            rawscale=max(D(1),abs(D(r['F'])),abs(2*D(r['G'])),abs(4*D(r['D'])))
            if abs(rawres-D(r['F_minus_2G_minus_4D']))>D('1e-8')*rawscale:fail('reported F-2G-4D mismatch')
        rows.append((z/q,dist/q,f/q,g/q,r['tau'],r['dead']))
    if seen!=set(rb):fail('missing globalrep row')
    M=D(len(rows));Z=sum(x[0] for x in rows)/M;DD=sum(x[1] for x in rows)/M;Q=DD/(D(N)*Z) if Z else Decimal('NaN')
    h=[x[1]-D(N)*Q*x[0] for x in rows];v=svar(h);se=(v/M).sqrt()/(D(N)*Z) if v is not None and Z else None
    print(json.dumps({'kind':'random_tau_aggregate','n':a.n,'K':K,'qden':qden,'observations':len(rows),'dead_observations':sum(x[5] for x in rows),'weighted_Z_mean':str(Z),'weighted_D_mean':str(DD),'ratio_Q':str(Q),'delta_SE_mixture':None if se is None else str(se),'weighted_F_mean':str(sum(x[2] for x in rows)/M),'weighted_G_mean':str(sum(x[3] for x in rows)/M),'tau_counts':{str(t):sum(x[4]==t for x in rows) for t in sorted({x[4] for x in rows})},'note':'Diagnostic delta SE only; it includes random-tau mixture and conditional SMC variance.'},indent=2))
if __name__=='__main__':main()
