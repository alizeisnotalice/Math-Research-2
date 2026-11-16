#!/usr/bin/env python3
"""Real-CLI fixtures for the log-scaled tau=2 distance aggregator."""
import copy, hashlib, json, math, pathlib, subprocess, sys, tempfile
from decimal import Decimal, getcontext

getcontext().prec=110
HERE=pathlib.Path(__file__).resolve().parent
AGG=HERE/'aggregate_logscaled_distance.py'
OLD_AGG=HERE/'frozen_r023_aggregate_tau2_distance.py'
SRC=HERE/'tau2_distance_terminal_r023.swift'
BASE=HERE/'n3-base.jsonl'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def near(a,b,eps=Decimal('1e-70')):
    assert abs(Decimal(str(a))-Decimal(str(b))) < eps, (a,b)

base=[json.loads(x) for x in BASE.read_text().splitlines() if x.strip()]
assert base[0]['kind']=='metadata' and base[-1]['kind']=='run_end'

def write_and_run(td, rows, epsilon='1/2', executable=AGG):
    raw=td/'fixture.jsonl'; raw.write_text('\n'.join(json.dumps(x) for x in rows)+'\n')
    m=rows[0]
    side=td/'sidecar.json'
    side.write_text(json.dumps({'raw_sha256':sha(raw),'source_sha256':sha(SRC),
        'expected_metadata':{k:m[k] for k in ('n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps')}}))
    p=subprocess.run([sys.executable,str(executable),str(raw),'--sidecar',str(side),'--source',str(SRC),'--epsilon',epsilon],text=True,capture_output=True)
    return p, raw, side

def consistent_reweight(r, z):
    """Rebuild every normalizer-bearing field for an exact fixture weight."""
    r['Z']=float(z)
    n=Decimal(base[0]['n']) ** 2
    u=Decimal(str(r['U_over_Z']))
    f0=Decimal(str(r['F_over_Z']))
    g0=Decimal(str(r['G_over_Z']))
    r['D']=float(z*n*u)
    r['F']=float(z*f0)
    r['G']=float(z*g0)
    r['F_minus_2G_minus_4D']=float(Decimal(str(r['F']))-2*Decimal(str(r['G']))-4*Decimal(str(r['D'])))
    for k in ('tau2_lift_complete_K2_Z','tau2_lift_complete_K2_D','tau2_lift_complete_K2_F','tau2_lift_complete_K2_G'):
        suffix=k.removeprefix('tau2_lift_complete_K2_')
        scale={'Z':'tau2_lift_scaled_mean_L','D':'tau2_lift_scaled_mean_L_distance2',
               'F':'tau2_lift_scaled_mean_L_F0','G':'tau2_lift_scaled_mean_L_G0'}[suffix]
        r[k]=float(z*Decimal(str(r[scale])))

def underflow_row(r, offset):
    x=copy.deepcopy(r)
    x['logZ']=str(Decimal(str(x['logZ']))-Decimal(offset))
    # This is precisely how the frozen Double producer reports values after a
    # true normalizer underflow; per-particle scale-free fields remain intact.
    for k in ('Z','D','F','G','F_minus_2G_minus_4D',
              'tau2_lift_complete_K2_Z','tau2_lift_complete_K2_D',
              'tau2_lift_complete_K2_F','tau2_lift_complete_K2_G'):
        x[k]=0
    return x

passed=[]
with tempfile.TemporaryDirectory(dir=HERE) as t:
    td=pathlib.Path(t)
    p,raw,side=write_and_run(td,base)
    assert p.returncode==0,p.stderr
    ordinary=json.loads(p.stdout); passed.append('ordinary frozen n3 input accepted')
    p,raw,side=write_and_run(td,base,executable=OLD_AGG)
    assert p.returncode==0,p.stderr
    oldordinary=json.loads(p.stdout)
    for k in ('Q','V','C','deltaSE_Q','distance_L_weighted_probabilities','distance_deltaSE','E2_tau','deltaSE_E2_tau','strict_tail_d_gt_epsilon_n','strict_tail_deltaSE'):
        if isinstance(ordinary[k],list):
            for a,b in zip(ordinary[k],oldordinary[k]): near(a,b)
        else: near(ordinary[k],oldordinary[k])
    passed.append('ordinary n3 new and frozen r023 aggregators agree')

    # The two real retained replicates are reweighted by 1:3.  A direct
    # hand calculation from the reported L-weighted numerators is compared
    # against the real aggregator CLI, including Q and the strict tail.
    rows=copy.deepcopy(base)
    rr=rows[1:-1]
    # Decimal strings make this a literal hand-built 1:3 normalizer fixture,
    # rather than merely a Double approximation to log(3).
    rr[0]['logZ']='0'
    rr[1]['logZ']=str(Decimal(3).ln())
    for r,z in zip(rr,(Decimal(1),Decimal(3))):
        consistent_reweight(r,z)
    p,raw,side=write_and_run(td,rows)
    assert p.returncode==0,p.stderr
    hand=json.loads(p.stdout)
    P=Decimal(rows[0]['particles'])
    lmean=lambda r: Decimal(sum(r['tau2_lift_L_groups']))/P
    dmean=lambda r: Decimal(sum(r['tau2_lift_L_distance2_groups']))/P
    den=sum(Decimal(z)*lmean(r) for z,r in zip((1,3),rr))
    qnum=sum(Decimal(z)*dmean(r) for z,r in zip((1,3),rr))
    near(hand['Q'],qnum/(Decimal(9)*den))
    inds=[d for d in range(5) if 2*d>3]
    hnum=sum(Decimal(z)*sum(Decimal(sum(g[d] for g in r['tau2_distance_L_hist_groups']))/P for d in inds) for z,r in zip((1,3),rr))
    near(hand['strict_tail_d_gt_epsilon_n'],hnum/den)
    passed.append('two-replicate hand weights 1:3 match CLI Q and strict tail')

    shifted=[rows[0]]+[underflow_row(r,1000) for r in rows[1:-1]]+[rows[-1]]
    p,raw,side=write_and_run(td,shifted)
    assert p.returncode==0,p.stderr
    shiftedout=json.loads(p.stdout)
    for k in ('Q','V','C','deltaSE_Q','distance_L_weighted_probabilities','distance_deltaSE','E2_tau','deltaSE_E2_tau','strict_tail_d_gt_epsilon_n','strict_tail_deltaSE','tau2_weight','liftZ_weight'):
        if isinstance(hand[k],list):
            for a,b in zip(hand[k],shiftedout[k]): near(a,b)
        elif isinstance(hand[k],dict):
            for q in hand[k]: near(hand[k][q],shiftedout[k][q],Decimal('1e-14'))
        else: near(hand[k],shiftedout[k])
    assert Decimal(shiftedout['log_lift_full_normalizer_mean']) < Decimal(-900)
    p,raw,side=write_and_run(td,shifted,executable=OLD_AGG)
    assert p.returncode==0,p.stderr
    oldshift=json.loads(p.stdout)
    for k in ('Q','V','C','deltaSE_Q','distance_L_weighted_probabilities','distance_deltaSE','E2_tau','deltaSE_E2_tau','strict_tail_d_gt_epsilon_n','strict_tail_deltaSE'):
        if isinstance(shiftedout[k],list):
            for a,b in zip(shiftedout[k],oldshift[k]): near(a,b)
        else: near(shiftedout[k],oldshift[k])
    passed.append('global logZ shift -1000 with true Double underflow is accepted by both aggregators and preserves ratios and SE')

    def rejected(name, mutate):
        x=copy.deepcopy(base); mutate(x); p,raw,side=write_and_run(td,x)
        assert p.returncode!=0,(name,p.stdout,p.stderr); passed.append(name+' rejected')
    rejected('NaN logZ',lambda x:x[1].__setitem__('logZ','NaN'))
    rejected('inconsistent nonzero Z/logZ',lambda x:x[1].__setitem__('logZ',float(x[1]['logZ'])+1))
    rejected('fake underflow',lambda x:x[1].__setitem__('Z',0))
    def badseed(x):
        x[1]['seed_hex']='0'
    rejected('bad deterministic seed',badseed)
    def badgroup(x):
        x[1]['tau2_lift_L_groups'][0]+=1
    rejected('bad integer group identity',badgroup)
    def bad_raw_scale_tuple(x):
        # Preserve F-2G=4D while making the entire reported root tuple
        # inconsistent with Z and its scale-free fields.
        for k in ('D','F','G'):
            x[1][k]*=2
    rejected('inconsistent raw root scale tuple',bad_raw_scale_tuple)
receipt={'passed':passed,'aggregator_sha256':sha(AGG),'frozen_r023_aggregator_sha256':sha(OLD_AGG),'source_sha256':sha(SRC),'test_sha256':sha(pathlib.Path(__file__))}
HERE.joinpath('test-logscaled-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
