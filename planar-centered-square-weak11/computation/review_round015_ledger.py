"""Read-only exact R15 fixed-position mass-transfer and ledger audit."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from collections import Counter
from itertools import product
from datetime import datetime,timezone
import argparse,hashlib,json,math,re,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def norm(x):return {'Px':str(F(x['Px'])),'Py':str(F(x['Py'])),'lam':str(F(x['lam'])),
                    'pts':[[str(F(v)) for v in q] for q in x['pts']]}
def apply(base,target,operations):
    n=len(base['pts']);orig=[F(p[2]) for p in target['pts']];weights=orig.copy()
    for op,coef in operations:
        ids=[op['quadrant']*n+j for j in op['members']];anchor=op['quadrant']*n+op['anchor']
        assert op['members']==[j for j,p in enumerate(base['pts']) if F(p[op['axis']])==F(op['coordinate'])]
        others=[j for j in range(n) if j not in op['members']]
        assert op['anchor']==max(others,key=lambda j:orig[op['quadrant']*n+j])
        G=sum(orig[j] for j in ids);delta=F(coef)*op['sign']*F(op['eps'])*min(G,orig[anchor])
        for j in ids:weights[j]+=delta*orig[j]/G
        weights[anchor]-=delta
    assert all(w>0 for w in weights) and sum(weights)==sum(orig)
    out=norm(target)
    for p,w in zip(out['pts'],weights):p[2]=str(w)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();start=time.perf_counter();t=a.topic;res=t/'computation/results/round015'
    ip=t/'computation/inputs/r015_final_0680.json';bp=t/'computation/inputs/r014_final_0935.json';r8p=t/'computation/inputs/r008_final_0611.json'
    d,target,base=map(read,(ip,bp,r8p));p,i,c=map(read,(res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json'))
    old=read(t/'computation/results/round014/primary_final.json');ledger=read(res/'search_ledger.json');n=len(base['pts'])
    assert n==51 and len(target['pts'])==len(d['pts'])==204
    assert sha(ip)==p['input_sha256']==i['input_sha256']==c['input_sha256'] and sha(bp)==c['baseline_r14_sha256'] and sha(r8p)==c['baseline_r8_sha256']
    assert all(F(u[k])==F(v[k]) for u,v in zip(target['pts'],d['pts']) for k in (0,1))
    assert all(F(d[k])==F(target[k]) for k in ('Px','Py','lam'))
    px,py,lam=map(F,(d['Px'],d['Py'],d['lam']));mass=sum(F(z[2]) for z in d['pts']);m8=sum(F(z[2]) for z in base['pts'])
    assert lam==1 and all(F(z[2])>0 for z in d['pts']) and mass==sum(F(z[2]) for z in target['pts'])==4*m8
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
    changed=[j for j,(u,v) in enumerate(zip(target['pts'],d['pts'])) if F(u[2])!=F(v[2])]
    assert c['changes_zero_based']=={'x':[],'y':[],'weight':changed} and changed==[1,32,33,103,134,135]
    for q in (0,2):
        for j in (32,33):assert F(d['pts'][q*n+j][2])==F(target['pts'][q*n+j][2])*F(10001,10000)
        assert F(target['pts'][q*n+1][2])-F(d['pts'][q*n+1][2])==sum(F(target['pts'][q*n+j][2]) for j in (32,33))/10000
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1
    assert tail==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    lo,up=F(p['inner_ratio']),F(p['outer_ratio']);gain=lo-F(old['outer_ratio'])
    assert lo==F(i['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(i['weak_ratio_upper']['fraction'])==F(c['upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(i['area_upper']['fraction'])==up*mass
    assert gain==F(c['gain_over_r14_upper'])>0 and F(p['tail_radius'])==F(i['radius_cap'])==10
    assert p['bits']==i['bits']==40 and p['code_sha256']==sha(t/'computation/periodic_halo.py') and i['script_sha256']==sha(t/'computation/independent_periodic_halo.py')
    assert sha(a.raw_search/'primary_old_R10.json')==sha(res/'primary_final.json')
    P=[F(base['Px']),F(base['Py'])];delta=F(0)
    for row,point in enumerate(target['pts']):
        quadrant,j=divmod(row,n);cx,cy=divmod(quadrant,2)
        for axis,copy in ((0,cx),(1,cy)):
            dd=(F(point[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(dd,2*P[axis]-dd))
    rho=m8/(P[0]*P[1]);original_weights=[F(z[2]) for z in target['pts']]
    all_hashes=set();actual=Counter();filemap={};batches=[];wall=Decimal(0);cpu=Decimal(0);chosen=None;raw_specs=None
    for index,row in enumerate(ledger):
        raw=a.raw_search/f'batch{index}';m=read(raw/'manifest.json');results=read(raw/'results.json');done=[json.loads(z) for z in (raw/'completed.jsonl').read_text().splitlines()]
        assert sha(raw/'manifest.json')==row['manifest_sha256'] and sha(raw/'results.json')==row['results_sha256']
        assert digest(m['tasks'])==m['task_list_sha256']==row['task_list_sha256']
        assert m['base_sha256']==sha(r8p) and m['target_sha256']==sha(bp) and F(m['delta'])==delta
        assert m['rejected']==[] and m['dominated_kernel_sha256']==row['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')
        assert len(results)==len(done)==len(m['tasks'])==row['points'] and len({z['id'] for z in results})==len(results)
        assert sorted(done,key=lambda z:z['id'])==sorted(results,key=lambda z:z['id'])
        tasks={z['id']:z for z in m['tasks']};rows={z['id']:z for z in results};filemap[row['results_sha256']]=(tasks,rows)
        alias_counter=Counter();coeff_set=set()
        for task in tasks.values():
            assert digest(task['input'])==task['input_sha256']
            for alias in task['aliases']:
                if index==0:
                    alias_counter[canon(alias)]+=1
                    operations=[] if alias.get('kind')=='baseline' else [(alias,F(1))]
                else:
                    cc=tuple(map(F,alias['coefficients']));coeff_set.add(cc);operations=list(zip(chosen,cc))
                assert apply(base,target,operations)==norm(task['input'])
        if index==0:
            expected=[{'kind':'baseline'}]
            for quadrant in range(4):
                for axis in (0,1):
                    groups={}
                    for j,z in enumerate(base['pts']):groups.setdefault(F(z[axis]),[]).append(j)
                    for value,ids in sorted(groups.items()):
                        others=[j for j in range(n) if j not in ids]
                        if not others:continue
                        anchor=max(others,key=lambda j:original_weights[quadrant*n+j])
                        for eps in ('1/1000','1/10000'):
                            for sign in (-1,1):expected.append({'quadrant':quadrant,'axis':axis,'members':ids,'coordinate':str(value),'anchor':anchor,'eps':eps,'sign':sign})
            assert Counter(map(canon,expected))==alias_counter;raw_specs=len(expected);assert raw_specs==529
            baseline=rows['0000'];rank=sorted((z for z in results if z['id']!='0000' and F(z['inner_ratio'])>F(baseline['outer_ratio'])),key=lambda z:F(z['inner_ratio']),reverse=True)
            chosen=[];seen=set()
            for z in rank:
                for alias in tasks[z['id']]['aliases']:
                    key=alias['quadrant'],tuple(alias['members']),alias['anchor']
                    if key not in seen:chosen.append(alias);seen.add(key);break
                if len(chosen)==2:break
            assert len(chosen)==2
        else:
            assert m['selection']['directions']==chosen
            assert coeff_set=={(F(x,4),F(y,4)) for x,y in product(range(-16,17),repeat=2)}
        new=sum(not z.get('reused',False) for z in results);reuse=len(results)-new;hashes={z['input_sha256'] for z in results}
        assert new==row['new_calls'] and reuse==row['reused'] and len(hashes-all_hashes)==row['new_distinct_inputs']
        for z in results:
            assert 'error' not in z and z['input_sha256']==tasks[z['id']]['input_sha256']
            td=tasks[z['id']]['input'];kappa=max(F(w)/F(base['pts'][j%n][2]) for j,(x,y,w) in enumerate(td['pts']))
            R=F(z['tail_radius']);bound=kappa*rho*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
            assert kappa==F(z['kappa']) and F(z['delta'])==delta and bound==F(z['tail_bound'])<1
            assert R==5 and F(z['inner_ratio'])<=F(z['outer_ratio'])
            if not z.get('reused',False):actual[z['input_sha256']]+=1
        win=max(results,key=lambda z:F(z['inner_ratio']));assert win['id']==row['winner']['id']
        if index==1:
            assert tasks[win['id']]['aliases']==[{'coefficients':['1','1']}]
            assert apply(base,target,list(zip(chosen,(F(1),F(1)))))==norm(d)
            assert F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        all_hashes|=hashes;wall+=Decimal(str(row['wall_seconds']));cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        batches.append({'batch':index,'rows':len(results),'actual':new,'reused':reuse,'tasks':tasks,'results':results})
    for batch in batches:
        for z in batch['results']:
            if not z.get('reused',False):continue
            tasks,rows=filemap[z['reused_results_sha256']];oldrow=rows[z['reused_from_id']]
            assert z['input_sha256']==oldrow['input_sha256'] and norm(batch['tasks'][z['id']]['input'])==norm(tasks[oldrow['id']]['input'])
            for key in ('kappa','delta','inner_ratio','outer_ratio'):assert F(z[key])==F(oldrow[key])
    total=sum(z['rows'] for z in batches);reuse=sum(z['reused'] for z in batches)
    assert (total,len(all_hashes),sum(actual.values()),len(actual),reuse)==(1586,1581,1581,1581,5)
    assert f'{wall:.6f}'=='690.211285' and f'{cpu:.6f}'=='6546.749079'
    frac=lambda x:[F(int(u),int(v)) for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',x)]
    proof=t/'proofs/round015-lower.md';report=t/'rounds/015.md';assert frac(proof.read_text())==[tail,lo,up,gain]
    r1=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/periodic-complete-halo.md');assert sha(r1)=='7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e'
    files=[ip,proof,report,res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','input_sha256':sha(ip),'baseline_R14_sha256':sha(bp),'positive_masses':204,'mass':str(mass),
         'changed_mass_indices_zero_based':changed,'coordinates_and_periods_unchanged':True,'chosen_directions':chosen,'winner_coefficients':['1','1'],
         'tail_radius':10,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R14_full_upper':str(gain),
         'ledger_counts':{'result_rows':total,'distinct_inputs':len(all_hashes),'actual':sum(actual.values()),'reused':reuse,'duplicate_actual':sum(actual.values())-len(actual)},
         'raw_first_batch_specs_before_deduplication':raw_specs,'first_batch_tasks_after_deduplication':497,
         'fixed_orig_mass_transfer_rule_independently_reconstructed':True,'saved_kappa_and_transport_bounds_exactly_recomputed':True,
         'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':i['total_seconds'],
         'geometry_or_search_rerun':False,'files_sha256':{str(p):sha(p) for p in files},'audit_script_sha256':sha(Path(__file__)),
         'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
