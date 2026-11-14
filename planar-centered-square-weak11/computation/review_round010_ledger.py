"""Read-only exact R10 supercell proof and deduplicated execution-ledger audit."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from collections import Counter
from datetime import datetime,timezone
import argparse,hashlib,json,math,platform,re,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def norm(x):return {'Px':str(F(x['Px'])),'Py':str(F(x['Py'])),'lam':str(F(x['lam'])),
                    'pts':[[str(F(v)) for v in q] for q in x['pts']]}
def make(base,alias):
    axis=alias['axis'];move=alias.get('move_axis',axis);P=[F(base['Px']),F(base['Py'])]
    periods=P.copy();periods[axis]*=2
    members=alias.get('members',[])
    shift=F(alias['shift']) if 'shift' in alias else alias.get('sign',0)*F(alias.get('eps','0'))
    if 'coordinate' in alias:
        assert members==[i for i,p in enumerate(base['pts']) if F(p[move])==F(alias['coordinate'])]
    pts=[]
    for copy in (0,1):
        for i,p in enumerate(base['pts']):
            q=list(map(F,p));q[axis]+=copy*P[axis]
            if copy==1 and i in members:q[move]+=shift*P[move]
            for k in (0,1):q[k]%=periods[k]
            pts.append(list(map(str,q)))
    return {'Px':str(periods[0]),'Py':str(periods[1]),'lam':str(F(base['lam'])),'pts':pts}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();start=time.perf_counter();t=a.topic;res=t/'computation/results/round010'
    ip=t/'computation/inputs/r010_final_0019.json';bp=t/'computation/inputs/r008_final_0611.json'
    d,b=read(ip),read(bp);p,i,c=map(read,(res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json'))
    old=read(t/'computation/results/round008/primary_final.json');ledger=read(res/'search_ledger.json')
    assert sha(ip)==p['input_sha256']==i['input_sha256']==c['input_sha256'] and sha(bp)==c['baseline_sha256']
    alias={'kind':'shift_second_copy','axis':1,'move_axis':1,'members':list(range(27,50,2)),'shift':'3/400000'}
    assert make(b,alias)==norm(d)
    mass=sum(F(q[2]) for q in d['pts']);oldmass=sum(F(q[2]) for q in b['pts'])
    px,py,lam=map(F,(d['Px'],d['Py'],d['lam']));assert lam==1
    assert len(d['pts'])==102 and mass==2*oldmass and all(F(q[2])>0 for q in d['pts'])
    assert px==F(b['Px']) and py==2*F(b['Py'])
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
    tail=mass/(px*py)*(1+px/16)*(1+py/16);q=math.ceil(F(8)/min(px,py))+1
    assert tail==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==4
    assert F(p['tail_radius'])==F(i['radius_cap'])==8 and p['bits']==i['bits']==40
    lo,up=F(p['inner_ratio']),F(p['outer_ratio']);margin=lo-F(old['outer_ratio'])
    assert lo==F(i['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(i['weak_ratio_upper']['fraction'])==F(c['upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lo*mass
    assert F(p['outer_area'])==F(i['area_upper']['fraction'])==up*mass
    assert margin==F(c['gain_over_baseline_upper'])>0
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py') and i['script_sha256']==sha(t/'computation/independent_periodic_halo.py')
    filemap={};batches=[];all_inputs=set();actual_hashes=Counter();values={};total=reused=early=0;wall=Decimal(0);cpu=Decimal(0)
    for row in ledger:
        raw=a.raw_search/row['path'];manifest=read(raw/'manifest.json');results=read(raw/'results.json')
        completed=[json.loads(z) for z in (raw/'completed.jsonl').read_text().splitlines()]
        assert sha(raw/'manifest.json')==row['manifest_sha256'] and sha(raw/'results.json')==row['results_sha256']
        assert digest(manifest['tasks'])==row['task_list_sha256']==manifest['task_list_sha256']
        assert manifest.get('source_input_sha256',manifest.get('input_sha256'))==sha(bp)
        assert len(results)==len(completed)==len(manifest['tasks'])==row['points']
        assert len({z['id'] for z in results})==len(results)
        assert sorted(completed,key=lambda z:z['id'])==sorted(results,key=lambda z:z['id'])
        tasks={z['id']:z for z in manifest['tasks']};byid={z['id']:z for z in results}
        filemap[row['results_sha256']]=(tasks,byid)
        for task in tasks.values():
            assert digest(task['input'])==task['input_sha256']
            assert all(make(b,alias)==norm(task['input']) for alias in task['aliases'])
        new=sum(not z.get('reused',False) for z in results);reuse=len(results)-new
        batch_hashes={z['input_sha256'] for z in results}
        assert new==row['new_worker_evaluations'] and reuse==row['reused']
        assert len(batch_hashes-all_inputs)==row['unique_new_across_earlier_batches']
        batch_early=0
        for z in results:
            assert 'error' not in z and z['input_sha256']==tasks[z['id']]['input_sha256'] and z['tail_radius']==8
            pair=F(z['inner_ratio']),F(z['outer_ratio']);assert pair[0]<=pair[1]
            if z['input_sha256'] in values:assert values[z['input_sha256']]==pair
            values[z['input_sha256']]=pair
            if not z.get('reused',False):
                td=tasks[z['id']]['input'];xx,yy,ll=map(F,(td['Px'],td['Py'],td['lam']));mm=sum(F(p[2]) for p in td['pts'])
                bounds=[mm/(xx*yy)*(1+xx/(2*r))*(1+yy/(2*r)) for r in range(5,9)]
                assert all(v>=ll for v in bounds[:3]) and bounds[-1]<ll
                batch_early+=3;actual_hashes[z['input_sha256']]+=1
        assert batch_early==row['tail_rejection_attempts']
        winner=max(results,key=lambda z:F(z['inner_ratio']));assert winner['id']==row['winner']['id']
        assert winner['input_sha256']==row['winner']['input_sha256'] and F(winner['inner_ratio'])==F(row['winner']['inner_ratio'])
        total+=len(results);reused+=reuse;early+=batch_early;all_inputs|=batch_hashes
        wall+=Decimal(str(row['wall_seconds']));cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        batches.append({'path':row['path'],'rows':len(results),'actual':new,'reused':reuse,'input_hashes':batch_hashes,'tasks':tasks,'results':results})
    reuse_links=0
    for batch in batches:
        for z in batch['results']:
            if not z.get('reused',False):continue
            current=z;visited=set()
            while current.get('reused',False):
                key=current['reused_results_sha256'],current['reused_from_id'];assert key not in visited;visited.add(key)
                tasks,rows=filemap[key[0]];source=rows[key[1]]
                assert source['input_sha256']==current['input_sha256']
                assert F(source['inner_ratio'])==F(current['inner_ratio']) and F(source['outer_ratio'])==F(current['outer_ratio'])
                assert norm(tasks[key[1]]['input'])==norm(batch['tasks'][z['id']]['input'])
                current=source
            reuse_links+=1
    actual=sum(actual_hashes.values());duplicates=actual-len(actual_hashes)
    assert (total,len(all_inputs),actual,reused,duplicates,early)==(358,296,326,32,30,978)
    assert len(actual_hashes)==len(all_inputs)==296 and reuse_links==32
    assert f'{wall:.6f}'=='279.958349' and f'{cpu:.6f}'=='2974.504144'
    assert all(c[k]==v for k,v in [('search_result_rows',total),('distinct_inputs',len(all_inputs)),('successful_worker_evaluations',actual),('reused_records',reused),('duplicated_worker_evaluations',duplicates),('early_tail_rejection_attempts',early)])
    frac=lambda x:[F(int(a.strip()),int(b.strip())) for a,b in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',x)]
    proof=t/'proofs/round010-lower.md';report=t/'rounds/010.md'
    assert frac(proof.read_text())==[tail,lo,up,margin]
    r1=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/periodic-complete-halo.md')
    assert sha(r1)=='7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e'
    files=[ip,proof,report,res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','input_sha256':sha(ip),'baseline_R8_sha256':sha(bp),'positive_masses':102,
         'mass':str(mass),'copy_factors':[1,2],'changed_copy_1_original_indices':list(range(27,50,2)),
         'exact_y_shift':str(3*F(b['Py'])/400000),'tail_radius':8,'tail_bound':str(tail),'q':q,
         'lower':str(lo),'upper':str(up),'lower_minus_R8_full_upper':str(margin),
         'ledger_counts':{'rows':total,'distinct_inputs':len(all_inputs),'successful_full_evaluations':actual,'reused':reused,'duplicate_full_evaluations':duplicates,'early_tail_only_rejections':early,'calculator_calls_including_early_rejections':actual+early},
         'batches':[{'path':z['path'],'rows':z['rows'],'actual':z['actual'],'reused':z['reused']} for z in batches],
         'reuse_chains_verified':reuse_links,'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),
         'geometry_or_search_rerun':False,'files_sha256':{str(p):sha(p) for p in files},'audit_script_sha256':sha(Path(__file__)),
         'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
