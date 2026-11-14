"""Read-only exact R6 input, saved-bound, proof and winning-ledger audit.
No geometry calculator, optimizer, or search is run.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from datetime import datetime,timezone
import argparse,copy,hashlib,json,math,platform,re,sys,time


def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(d):return json.dumps(d,sort_keys=True,separators=(',',':')).encode()
def digest(d):return hashlib.sha256(canonical(d)).hexdigest()
def normalize(d):return {'lam':str(F(d['lam'])),'Px':str(F(d['Px'])),'Py':str(F(d['Py'])),
                         'pts':[[str(F(v)) for v in p] for p in d['pts']]}


def apply_alias(base,alias):
    out=normalize(base);pts=[[F(v) for v in p] for p in out['pts']]
    g=alias['group'];axis=g['axis'];members=g['members'];coordinate=F(g['coordinate'])
    assert members==[i for i,p in enumerate(pts) if p[axis]==coordinate]
    oldmass=sum(p[2] for p in pts);eps=F(alias['eps'])*alias['sign']
    if alias['kind']=='coordinate':
        period=F(out['Px' if axis==0 else 'Py'])
        for i in members:pts[i][axis]=(pts[i][axis]+period*eps)%period
        assert alias['anchor'] is None
    else:
        assert alias['kind']=='mass'
        others=[i for i in range(len(pts)) if i not in members]
        anchor=max(others,key=lambda i:pts[i][2]);assert anchor==alias['anchor']
        groupmass=sum(pts[i][2] for i in members)
        delta=eps*min(groupmass,pts[anchor][2])
        for i in members:pts[i][2]+=delta*pts[i][2]/groupmass
        pts[anchor][2]-=delta
    assert all(p[2]>0 for p in pts) and sum(p[2] for p in pts)==oldmass
    out['pts']=[[str(v) for v in p] for p in pts]
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.perf_counter();topic=args.topic
    inp=topic/'computation/inputs/r006_final_0015.json';d=read(inp)
    res=topic/'computation/results/round006';p=read(res/'primary_final.json');i=read(res/'independent_final.json')
    check=read(res/'integration_checks.json');arithmetic=read(res/'arithmetic_check.json')
    baseline=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/computation/inputs/r004_final_0504.json')
    b=read(baseline);bresult=read(Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round004/independent_final.json'))
    inputsha=sha(inp);assert inputsha==p['input_sha256']==i['input_sha256']==check['input_sha256']
    assert sha(baseline)==check['baseline_input_sha256']
    assert len(d['pts'])==len(b['pts'])==51 and all(F(q[2])>0 for q in d['pts'])
    px,py,lam=F(d['Px']),F(d['Py']),F(d['lam'])
    assert all(F(d[k])==F(b[k]) for k in ('Px','Py','lam')) and lam==1
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,m in d['pts'])
    mass=sum((F(q[2]) for q in d['pts']),F(0));assert mass==sum((F(q[2]) for q in b['pts']),F(0))
    tail=mass/(px*py)*(1+px/10)*(1+py/10);q=math.ceil(F(5)/min(px,py))+1
    assert tail<1 and q==3 and tail==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])==F(check['tail_bound'])
    assert mass==F(p['mass'])==F(i['base_mass'])==F(check['mass'])
    lower,upper=F(p['inner_ratio']),F(p['outer_ratio'])
    assert lower==F(i['weak_ratio_lower']['fraction'])==F(check['lower'])==F(arithmetic['lower'])
    assert upper==F(i['weak_ratio_upper']['fraction'])==F(check['upper'])==F(arithmetic['seed_upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lower*mass
    assert F(p['outer_area'])==F(i['area_upper']['fraction'])==upper*mass
    margin=lower-F(bresult['weak_ratio_upper']['fraction'])
    assert margin==F(check['difference_over_baseline_outer'])==F(arithmetic['gain_over_previous_seed_upper'])>0
    assert p['bits']==i['bits']==check['bits']==40 and F(check['epsilon'])==F(1,2**41)
    assert i['script_sha256']==sha(topic/'computation/independent_periodic_halo.py')
    assert p['code_sha256']==sha(topic/'computation/periodic_halo.py')
    assert all(i['primary_comparison']['checks'].values())
    changes=[];indexsets={k:[] for k in ('x','y','weight')}
    for j,(a,c) in enumerate(zip(b['pts'],d['pts'])):
        for k,name in enumerate(('x','y','weight')):
            if F(a[k])!=F(c[k]):
                indexsets[name].append(j)
                changes.append({'zero_based_atom':j,'field':name,'before':str(F(a[k])),
                                'after':str(F(c[k])),'delta':str(F(c[k])-F(a[k]))})
    assert indexsets==check['coordinate_changes_indices_zero_based']
    assert indexsets=={'x':[2,3,4,5],'y':list(range(27,50,2)),'weight':list(range(6))}
    ledger=read(res/'search_ledger.json');assert len(ledger)==8
    current=normalize(b);previous_outer=F(bresult['weak_ratio_upper']['fraction']);chain=[]
    wall=Decimal(0);cpu=Decimal(0);total=0
    for epoch,row in enumerate(ledger):
        summary=read(res/f'epoch{epoch}_summary.json')
        raw=args.raw_search/f'epoch{epoch}'
        manifest=read(raw/'manifest.json');results=read(raw/'results.json')
        completions=[json.loads(line) for line in (raw/'completed.jsonl').read_text().splitlines()]
        assert row['epoch']==summary['epoch']==epoch
        assert row['points']==row['completed']==summary['points']==summary['completed']==513
        assert row['errors']==summary['errors']==0 and row['workers']==12 and row['groups']==33
        assert row['bits']==40
        assert sha(raw/'manifest.json')==row['manifest_sha256']
        assert sha(raw/'results.json')==row['results_sha256']==summary['results_sha256']
        assert digest(manifest['tasks'])==row['task_list_sha256']==summary['manifest_task_list_sha256']
        assert digest(current)==row['baseline_canonical_sha256']==manifest['baseline_canonical_sha256']
        assert F(row['baseline_outer_ratio'])==previous_outer
        assert len(manifest['tasks'])==len(results)==len(completions)==513
        assert len({r['id'] for r in results})==len({r['id'] for r in completions})==513
        assert sorted(completions,key=lambda r:r['id'])==sorted(results,key=lambda r:r['id'])
        tasks={t['id']:t for t in manifest['tasks']}
        assert all(digest(t['input'])==t['input_sha256'] for t in tasks.values())
        assert all('error' not in r and r['input_sha256']==tasks[r['id']]['input_sha256']
                   and F(r['inner_ratio'])<=F(r['outer_ratio']) and r['tail_radius']==5 for r in results)
        winner=max(results,key=lambda r:F(r['inner_ratio']))
        top=summary['top'][0]
        assert winner['id']==row['winner_id']==top['id']
        assert F(winner['inner_ratio'])==F(row['winner_inner_ratio'])==F(top['inner_ratio'])>previous_outer
        candidates=[apply_alias(current,alias) for alias in row['winner_aliases']]
        assert all(c==candidates[0] for c in candidates)
        current=candidates[0]
        assert digest(current)==row['winner_input_canonical_sha256']==top['input_sha256']
        assert normalize(tasks[winner['id']]['input'])==current
        winnerfile=raw/f"candidate_{winner['id']}.json"
        assert normalize(read(winnerfile))==current
        if epoch<7:assert sha(winnerfile)==ledger[epoch+1]['baseline_file_sha256']
        previous_outer=F(winner['outer_ratio'])
        assert previous_outer==F(top['outer_ratio'])
        wall+=Decimal(str(row['wall_seconds']))
        cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        total+=len(results)
        chain.append({'epoch':epoch,'winner_id':winner['id'],'canonical_input_sha256':digest(current),
                      'strict_gain':str(F(winner['inner_ratio'])-F(row['baseline_outer_ratio'])),
                      'completed_records':len(completions),'completion_sha256':sha(raw/'completed.jsonl'),
                      'aliases_replayed':len(candidates),'final_aliases':row['winner_aliases']})
    assert current==normalize(d) and previous_outer==upper and F(ledger[-1]['winner_inner_ratio'])==lower
    assert total==4104 and f'{wall:.6f}'=='666.976484' and f'{cpu:.6f}'=='7095.199049'
    proof=topic/'proofs/round006-lower.md';report=topic/'rounds/006.md'
    fractions=lambda s:[F(int(a.strip()),int(b.strip())) for a,b in re.findall(r'\\frac\{([^{}]+)\}\s*\{([^{}]+)\}',s)]
    assert fractions(proof.read_text())==[mass,tail,lower,upper,margin]
    assert fractions(report.read_text())==[lower]
    r1=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/periodic-complete-halo.md')
    assert sha(r1)=='7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e'
    files=[inp,proof,report,res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','scope':'exact saved arithmetic, actual input and replayed ledger; no geometry/search rerun',
         'verified_utc':datetime.now(timezone.utc).isoformat(),'python':sys.version,'platform':platform.platform(),
         'input_sha256':inputsha,'positive_masses':51,'mass':str(mass),'tail_bound':str(tail),'q':q,
         'lower':str(lower),'upper':str(upper),'lower_minus_round004_full_upper':str(margin),
         'change_indices_zero_based':indexsets,'exact_changes':changes,'winning_transform_replay':chain,
         'search_completed_records_verified':total,'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),
         'files_sha256':{str(p):sha(p) for p in files},
         'audit_script_sha256':sha(Path(__file__)),'elapsed_seconds':time.perf_counter()-start,
         'command':[sys.executable,*sys.argv]}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('exact_changes','winning_transform_replay','files_sha256')},indent=2))

if __name__=='__main__':main()
