"""Audit old outputs, dense rank shortcut, sparse n11 and production refusals."""
import argparse
import json
from math import factorial, prod
from pathlib import Path
import subprocess

from count_protected import HERE, MIB, command, count, digest, frequencies, plan, verified_primes


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('old_histograms',type=Path)
    parser.add_argument('old_pratt',type=Path)
    a = parser.parse_args()
    primes,nodes=verified_primes(HERE/'pratt-61bit.json')
    new=read(HERE/'pratt-61bit.json')
    old=read(a.old_pratt)
    assert new['roots'][:4]==old['roots'] and all(new['nodes'][p]==v for p,v in old['nodes'].items())
    assert nodes==50 and primes[-1]==2**61-259
    assert factorial(65)<prod(primes)<factorial(66)
    replay_file=HERE/'small-replay'/'summary.json'
    replay=read(replay_file)
    ref=read(a.old_histograms)
    assert replay['completed_jobs']==20
    names=('counts','unsafe_max','unsafe_second_after_safe_nonadjacent_max')
    for h,rows in replay['region_histograms'].items():
        for field,actual in zip(names,rows):
            expected=ref['regions'][h][field]
            m=max(len(actual),len(expected))
            assert actual+['0']*(m-len(actual)) == expected+['0']*(m-len(expected)), (h,field)
    for row in replay['upper_counts']:
        old=next(s for s in ref['upper_bounds'] if s['profile']==row['cut'])
        assert row['A2_upper_count']==old['protected_top_two_refined_upper_A2']
        assert row['global_second_excluded_count']==old['unprotected_nonadjacent_global_second_removed_words']
    # The DP's direct-rank mode checks every removal shortcut and the
    # enumeration rank, including >1 million protected n9 states.
    audit_out=[]
    for h in ('112','2322','444445555','555554444'):
        cmd=command(h,4,32)+['--verify-ranks']
        p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=30)
        assert not p.stderr
        r=json.loads(p.stdout)
        assert r['states_by_layer']==r['processed_by_layer']
        if h in replay['region_histograms']:
            old_r=read(HERE/'small-replay'/f'{h}-p4.json')
            for k in ('peak_frequency_residues','unprotected_max_frequency_residues','bad_second_frequency_residues'):
                assert r[k]==old_r[k]
        audit_out.append({'profile':h,'states':r['independent_sets'],
                          'verified_deleted_predecessor_ranks':sum(k*s for k,s in enumerate(r['states_by_layer']))})
    # A path with one extra vertex supports a sole protected point of degree 2.
    h='11111111112'
    sparse=HERE/'sparse-n11-gate';sparse.mkdir(exist_ok=False)
    plan_sparse=plan(h,1,primes)
    residues=[count(h,ix,p,sparse,plan_sparse,1,10) for ix,p in enumerate(primes)]
    actual=frequencies(h,{key:r for key,r,_ in residues},primes)
    m=sum(map(int,h))
    peak=factorial(m)//3
    unsafe=(m-3)*factorial(m-1)//3
    second=(m-3)*factorial(m-2)
    assert actual==([factorial(m)-peak,peak], [(m-1)*factorial(m-1)-unsafe,unsafe], [0,second])
    near=('55555566666','66666655555')
    mirror_plans=[plan(q,6656,primes) for q in near]
    assert [x['max_table_bytes'] for x in mirror_plans]==[4389077248,6589688840]
    assert all(r['fits_budget'] and r['five_moduli_cover_factorial'] for r in mirror_plans)
    reject=subprocess.run(command(near[1],0,6284),capture_output=True,text=True,timeout=5)
    assert reject.returncode==2 and not reject.stdout and 'two-table budget refused' in reject.stderr
    no_crt=HERE/'crt-refusal-gate'
    reject=subprocess.run(['python3','-B',str(HERE/'count_protected.py'),str(no_crt),'55555555555','--execute'],
                          capture_output=True,text=True,timeout=5)
    assert reject.returncode==2 and 'five moduli do not cover' in reject.stderr
    assert set(z.name for z in no_crt.iterdir())=={'plan.json'}
    report={'counter_sha256':digest(HERE/'protected_mod.cpp'),'binary_sha256':digest(HERE/'protected_mod'),
            'old_pratt_sha256':digest(a.old_pratt),'old_histograms_sha256':digest(a.old_histograms),
            'new_pratt_sha256':digest(HERE/'pratt-61bit.json'),'prior_histogram_replay':sorted(replay['region_histograms']),
            'replay_jobs_waited':replay['completed_jobs'],'direct_rank_gates':audit_out,'sparse_n11_histograms':actual,
            'sparse_receipts':[r for _,_,r in residues],'lower_budget_refused':True,'straight_exact_66_refused':True,
            'balance_plans':mirror_plans,'scope':'Exact old-ledger and rank gates, no full n11 locally counted'}
    (HERE/'gate-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('exact old ledger, direct rank, sparse and allocation-refusal gates passed', audit_out)

if __name__=='__main__':
    main()
