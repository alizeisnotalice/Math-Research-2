"""Audit all four scored stages and transport whole-box bounds to the next target."""
from argparse import ArgumentParser
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np

from audit_tomar_outer_083730 import AREA, sha


def pick(path,*names):
    with np.load(path,allow_pickle=False) as z:
        return [z[n] for n in names]


def native_ids(path):
    info=json.loads(path.read_bytes())
    groups=info.get('retained',info)
    return {(int(p),i) for p,r in groups.items() for i in r}


def main():
    p=ArgumentParser(description=__doc__)
    for n in ('work','first_seeds','first_seed_audit','future_seeds','future_seed_audit'):
        p.add_argument(n,type=Path)
    p.add_argument('--old-inner',type=Path)
    p.add_argument('--old-outer',type=Path)
    p.add_argument('--next-emit',type=Path)
    a=p.parse_args()
    run_raw=(a.work/'stage-ledger.json').read_bytes()
    run=json.loads(run_raw)
    initial=json.loads(a.first_seed_audit.read_bytes())
    next_=json.loads(a.future_seed_audit.read_bytes())
    assert run['seed_receipt_sha']==sha(a.first_seed_audit.read_bytes())
    assert run['target']==initial['target'] and run['radial_target']==initial['radial_target']
    assert run['integer_goal']==initial['target_integer']
    assert run['complete'] and run['parts']==1024 and len(run['audits'])==4
    assert initial['radii']==next_['radii'] and initial['target_integer']<next_['target_integer']
    assert initial['provenance']==next_['provenance']
    goal=run['integer_goal'];higher=next_['target_integer']
    assert F(goal,AREA)>F(initial['target']) and F(higher,AREA)>F(next_['target'])
    parts=[{r[0]:r[1:] for r in item['rows']} for item in (initial,next_)]
    assert set(parts[0])==set(parts[1]) and len(parts[0])==1024
    scores_by_stage=tuple({p:h for p,h in k['work_files']} for k in run['audits'])
    assert all(set(z)==set(parts[0]) for z in scores_by_stage)
    accepted_next=0
    next_count=0
    retained={}
    archive=[]
    totals=[dict(tested=0,accepted=0,remaining=0,at_next_target=0,capped_by_radial=0,calls=0)
            for _ in range(4)]
    if a.next_emit:
        a.next_emit.mkdir(parents=True,exist_ok=True)
    for part,row in sorted(parts[0].items()):
        name=f'part{part:04}.npz'
        first=a.first_seeds/name
        last_sha=sha(first.read_bytes())
        active,old_score=pick(first,'pos','score')
        assert [sha(active.tobytes()),sha(old_score.tobytes())]==row[-2:]
        future,future_score=pick(a.future_seeds/name,'pos','score')
        assert [sha(future.tobytes()),sha(future_score.tobytes())]==parts[1][part][-2:]
        assert np.isin(active,future,assume_unique=True).all()
        to_skip=[]
        for j in range(4):
            path=a.work/f'stage{j+1}'/name
            payload_sha=sha(path.read_bytes())
            record=json.loads(path.with_suffix('.json').read_bytes())
            config=record['config']
            assert config['part']==part and config['layer']==j+1 and config['goal']==goal
            assert config['radial_goal']==run['radial_target'] and config['kernel_sha256']==run['source_sha']
            assert config['source_sha256']==last_sha and record['data_sha256']==payload_sha==scores_by_stage[j][part]
            last_sha=payload_sha
            passed,remaining,full,clipped=pick(path,'accepted','remaining','full_floor','clipped_floor')
            assert active.dtype==passed.dtype==remaining.dtype==np.int32
            assert np.array_equal(np.sort(np.r_[passed,remaining]),active)
            assert full.dtype==clipped.dtype==np.int64 and full.shape==clipped.shape==passed.shape
            assert np.all(full>=goal) and np.all((clipped==0)|(clipped>=full))
            limited=(clipped==0)|(clipped>full)
            assert np.all(np.isin(full[limited],[x['integer_area_lower'] for x in initial['radii']]))
            strong=passed[full>=higher]
            c=dict(tested=len(active),accepted=len(passed),remaining=len(remaining),calls=record['counts'][0],
                   at_next_target=len(strong),capped_by_radial=int(((clipped!=0)&(full<clipped)).sum()))
            assert all(record[k]==v for k,v in c.items() if k!='calls')
            assert record['counts'][2]==len(passed) and record['counts'][3]==len(remaining)
            for k,v in c.items():
                totals[j][k]+=v
            to_skip.extend(strong.tolist())
            active=remaining
        if len(active):
            retained[str(part)]=active.tolist()
        assert np.isin(np.array(to_skip,dtype=np.int32),future,assume_unique=True).all()
        keep=~np.isin(future,np.array(to_skip,dtype=np.int32),assume_unique=True)
        accepted_next+=len(to_skip)
        next_count+=int(keep.sum())
        archive.append([part,len(future),len(to_skip),int(keep.sum()),sha(future[keep].tobytes())])
        if a.next_emit:
            np.savez_compressed(a.next_emit/name,pos=future[keep],score=future_score[keep])

    for old,new in zip(run['audits'],totals,strict=True):
        assert all(old[k]==v for k,v in new.items())
    assert retained==json.loads((a.work/'retained.json').read_bytes())
    assert sum(map(len,retained.values()))==totals[-1]['remaining']
    assert next_count+accepted_next==next_['totals']['remaining']
    result=dict(stage_sha256=sha(run_raw),first_seed_sha256=sha(a.first_seed_audit.read_bytes()),
           next_seed_sha256=sha(a.future_seed_audit.read_bytes()),stages=totals,
           target=run['target'],next_target=next_['target'],
           dropped_by_scores_at_next=accepted_next,still_to_check_at_next=next_count,
           next_target_rows=archive,retained=retained,script_sha256=sha(Path(__file__).read_bytes()))
    if a.old_inner and a.old_outer:
        A=native_ids(a.old_inner); B=native_ids(a.old_outer)
        assert len(A)==501 and len(B)==40 and not (A&B)
        new={(int(p),i) for p,T in retained.items() for i in T}
        extra=new-(A|B)
        D={}
        for p,i in sorted(extra):
            D.setdefault(str(p),[]).append(i)
        result['previous_541']={'intersection_count':len(new&(A|B)),'new_count':len(extra),'new_leaves':D}
    print(json.dumps(result,separators=(',',':')))


if __name__=='__main__':
    main()
