"""Run two frozen Tomar guards on a preregistered family of entire original leaves."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as F
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from flint import arb

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_paired_affine_cap import digest,compile_drivers
from replay_paired_affine_tree import paired_driver


def freeze_cases(plan, prior_rows):
    def select(part, low=0, high=10**12, point=True):
        records=[r for r in prior_rows if r[0]==part and low<=r[3]<=high
                 and (r[5]>0 if point else r[5]==0)]
        return min(records,key=lambda r:(-r[5],-r[3],r[1]))
    expected=[
        ('highest_P30_among_176_between_250000_and_1200000_prior_nodes',select(176,250000,1200000)),
        ('mandated_180_24225',next(r for r in prior_rows if r[:2]==[180,24225])),
        ('highest_P30_in_181',select(181)),
        ('highest_P30_in_1205',select(1205)),
        ('highest_prior_node_count_WO_P30_in_1205',select(1205,point=False)),
    ]
    assert [(r['rule'],r['prior']) for r in plan['cohort']]==expected
    return [r for _,r in expected]


def summary(results):
    fields={}
    for kind in ('original','paired'):
        done=[(r['key'],r[kind]) for r in results]
        closed=[(key,stats[1]) for key,stats in done if stats[0]]
        fields[kind]=dict(
            closed=len(closed),
            worst_observed=max(stats[1] for _,stats in done),
            worst_closed=max((count for _,count in closed),default=None),
            closed_nodes=sum(count for _,count in closed),
            visit_sum=sum(stats[1] for _,stats in done),
            unresolved=[key for key,stats in done if not stats[0]]
        )
    fields['paired_accepted_total']=sum(r['paired_counters'][0] for r in results)
    closed=[r for r in results if r['paired'][0] and r['original'][0]]
    fields['both_closed_count']=len(closed)
    fields['both_closed_old_nodes']=sum(r['original'][1] for r in closed)
    fields['both_closed_new_nodes']=sum(r['paired'][1] for r in closed)
    return fields


def run(a):
    sys.path.insert(0,str(a.comp.resolve()))
    from certify_tomar_six_high import FINE,Q,SourceTree,build_kernel,load_table,ceil
    from replay_finite_point_083732 import PointProof
    from certify_nine import ar
    plan=json.loads(a.plan.read_bytes())
    assert F(plan['target'])==F('.83744') and F(plan['prior_target'])==F('.83742')
    assert plan['node_limit']==1500000 and plan['depth_limit']==135
    ledger=a.comp/'finite-point-lift-083742.json'
    assert digest(ledger)==plan['prior_sha256']=='0e84a6fe33567f0784f08b88db7cae6dc79ef6807891bc175d510900139c8cec'
    data=json.loads(ledger.read_bytes())
    cases=freeze_cases(plan,data['rows'])
    for name,key in (('tomar-six-high-input.json','base_sha256'),
                     ('finite-point-leaf-input.json','point_sha256'),
                     ('finite-point-global-input.json','bank_sha256')):
        assert digest(a.comp/name)==data[key]
    base=json.loads((a.comp/'tomar-six-high-input.json').read_bytes())
    cloud=json.loads((a.comp/'finite-point-leaf-input.json').read_bytes())
    bank=json.loads((a.comp/'finite-point-global-input.json').read_bytes())
    proof=PointProof(cloud,bank['cells'],100000000008865717121)
    pi=arb.pi()
    radius=ar(F(707,1000))
    assert (radius*radius<ar(F(1,2)) and
            (pi-(1/(2*radius)).acos())/4+(radius*radius-ar(F(1,4))).sqrt()/2>ar(F('.838')))

    a.work.mkdir(parents=True,exist_ok=True)
    _,original,fan_sha=compile_drivers(a.comp,a.work,a.cache,base,data,build_kernel)
    paired,source_sha=paired_driver(a.comp,a.work,a.cache)
    tree=SourceTree(a.tree)
    table={n:load_table(getattr(a,f'table{n}'),tree,base['tables'][str(n)])
           for n in sorted({r[2] for r in cases})}
    origin=np.array([int(F(x)*FINE) for x in tree.root['hi']],np.int64)
    goal=ceil(F(plan['target'])*2*Q*Q)
    assert goal==120687823102244742 and F(goal,2*Q*Q)>F(plan['target'])
    results={}
    timings={}

    def trial(key,mode,bounds,support,tbl):
        start=monotonic()
        print(f'{key} {mode} starts',file=sys.stderr,flush=True)
        kernel=original if mode=='original' else paired
        lo,hi=bounds
        status=np.zeros(34,np.int64)
        assert not kernel(lo,hi,origin,lo,hi,*tbl,goal,plan['depth_limit'],plan['node_limit'],
                          *support,len(support[1]),len(support[3])-1,status)
        stats=status[:7].tolist()
        assert 0<stats[1]<=plan['node_limit']
        if stats[0]:
            assert stats[5]==0 and stats[1]==2*(stats[2]+stats[3])-1
        else:
            assert stats[5] in (1,2,3,4)
            assert stats[5]!=1 or stats[1]==plan['node_limit']
        row={mode:stats}
        if mode=='paired':
            row['paired_counters']=status[7:12].tolist()
            if not stats[0]:
                row['first_unresolved_box']=[status[12:23].tolist(),status[23:34].tolist()]
        sec=monotonic()-start
        print(f'{key} {mode} ends {sec:.3f}s {stats}, paired extra {row.get("paired_counters")}',
              file=sys.stderr,flush=True)
        return key,mode,row,sec

    prepared={}
    for part,pos,level,*_ in cases:
        key=(part,pos)
        box=tree.leaf(*key)
        bounds=[np.array([int(x*FINE) for x in xrow],np.int64) for xrow in zip(*box,strict=True)]
        support=proof.presolve(box)
        assert len(support[3])-1==384
        prepared[key]=(bounds,support,table[level])
        results[key]=dict(key=list(key),material_level=level,root=[x.tolist() for x in bounds])
        timings[key]={}
    # Start paired and original work independently; every job has its own full bounds.
    order=sorted(prepared,key=lambda key:(key== (180,24225),key))
    with ThreadPoolExecutor(a.workers) as pool:
        jobs=[pool.submit(trial,key,mode,*prepared[key])
              for key in order for mode in ('paired','original')]
        for job in as_completed(jobs):
            key,mode,row,sec=job.result()
            results[key].update(row)
            timings[key][mode]=sec
    if a.times:
        a.times.write_text(json.dumps([{'key':list(key),**timings[key]} for key in sorted(timings)],
                                      separators=(',',':'))+'\n')
    frozen=[]
    for key in sorted(results):
        r=results[key]
        z={n:r[n] for n in ('key','material_level','root','original','paired','paired_counters')}
        if 'first_unresolved_box' in r:
            z['first_unresolved_box']=r['first_unresolved_box']
        frozen.append(z)
    return dict(target=plan['target'],target_integer=goal,round182_sha256=digest(ledger),
                plan_sha256=digest(a.plan),driver_sha256=digest(Path(__file__)),
                paired_source_sha256=source_sha,fan_sha256=fan_sha,
                depth_limit=plan['depth_limit'],node_limit=plan['node_limit'],
                results=frozen,summary=summary(frozen))


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for name in ('comp','tree','table32','table64','work','cache','plan'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--times',type=Path)
    a=parser.parse_args()
    assert 1<=a.workers<=4
    print(json.dumps(run(a),separators=(',',':')))
