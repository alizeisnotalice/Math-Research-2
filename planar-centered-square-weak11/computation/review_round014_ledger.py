"""Read-only exact R14 mode/grid, saved-certificate and reuse-ledger audit."""
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
def make(base,coeff):
    A,B,C=map(F,coeff);px,py=F(base['Px']),F(base['Py']);pts=[];delta=F(0);group=set(range(27,50,2))
    for x in (0,1):
        for y in (0,1):
            for i,p in enumerate(base['pts']):
                xx,yy,w=map(F,p);dy=F(0)
                if y==1 and i in group:dy+=A*3*py/400000
                if x==y==1 and i==50:dy+=B*py/400000
                if y==0 and i in group:dy+=C*(2*x-1)*py/100000
                delta=max(delta,abs(dy));pts.append([str((xx+x*px)%(2*px)),str((yy+y*py+dy)%(2*py)),str(w)])
    return {'Px':str(2*px),'Py':str(2*py),'lam':str(F(base['lam'])),'pts':pts},delta

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();start=time.perf_counter();t=a.topic;res=t/'computation/results/round014'
    ip=t/'computation/inputs/r014_final_0935.json';bp=t/'computation/inputs/r012_final_0017.json';r8p=t/'computation/inputs/r008_final_0611.json'
    d,b,r8=map(read,(ip,bp,r8p));p,i,c=map(read,(res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json'))
    old=read(t/'computation/results/round012/primary_final.json');ledger=read(res/'search_ledger.json')
    assert sha(ip)==p['input_sha256']==i['input_sha256']==c['input_sha256']
    assert sha(bp)==c['baseline_r12_sha256'] and sha(r8p)==c['baseline_r8_sha256']
    coeff=[F(19,16),F(27,64),F(0)];expected,delta=make(r8,coeff)
    assert norm(d)==expected and list(map(F,c['coefficients']))==coeff
    mass=sum(F(z[2]) for z in d['pts']);m8=sum(F(z[2]) for z in r8['pts']);px,py,lam=map(F,(d['Px'],d['Py'],d['lam']))
    assert len(d['pts'])==204 and all(F(z[2])>0 for z in d['pts']) and mass==4*m8
    assert (px,py,lam)==(2*F(r8['Px']),2*F(r8['Py']),F(1))
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
    changes={k:[] for k in ('x','y','weight')}
    for j,(u,v) in enumerate(zip(b['pts'],d['pts'])):
        for axis,name in enumerate(('x','y','weight')):
            if F(u[axis])!=F(v[axis]):changes[name].append(j)
    assert changes==c['changes_zero_based']
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1
    assert tail==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    assert F(p['tail_radius'])==F(i['radius_cap'])==10 and p['bits']==i['bits']==40
    lo,up=F(p['inner_ratio']),F(p['outer_ratio']);gain=lo-F(old['outer_ratio'])
    assert lo==F(i['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(i['weak_ratio_upper']['fraction'])==F(c['upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(i['area_upper']['fraction'])==up*mass
    assert gain==F(c['gain_over_r12_upper'])>0
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py') and i['script_sha256']==sha(t/'computation/independent_periodic_halo.py')
    assert sha(a.raw_search/'primary_old_R10.json')==sha(res/'primary_final.json')
    raw_ind=a.raw_search.parent/'independent_final.json';assert sha(raw_ind)==sha(res/'independent_final.json')
    comparison=read(res/'comparison.json');assert comparison['primary_file_sha256']==sha(res/'primary_final.json') and comparison['independent_file_sha256']==sha(res/'independent_final.json')
    all_hashes=set();actual=Counter();batches=[];filemap={};wall=Decimal(0);cpu=Decimal(0);all_coeff=[]
    rho=m8/(F(r8['Px'])*F(r8['Py']));bycoeff={};coarse_win=None
    for index,row in enumerate(ledger):
        raw=a.raw_search/f'batch{index}';m=read(raw/'manifest.json');results=read(raw/'results.json')
        done=[json.loads(z) for z in (raw/'completed.jsonl').read_text().splitlines()]
        assert sha(raw/'manifest.json')==row['manifest_sha256'] and sha(raw/'results.json')==row['results_sha256']
        assert digest(m['tasks'])==m['task_list_sha256']==row['task_list_sha256']
        assert m['input_sha256']==sha(r8p) and m['copy_factors']==[2,2] and m['tail_radius']==5
        assert len(results)==len(done)==len(m['tasks'])==row['points'] and len({z['id'] for z in results})==len(results)
        assert sorted(done,key=lambda z:z['id'])==sorted(results,key=lambda z:z['id'])
        tasks={z['id']:z for z in m['tasks']};rows={z['id']:z for z in results};filemap[row['results_sha256']]=(tasks,rows)
        coeff_set=set()
        for task in tasks.values():
            assert digest(task['input'])==task['input_sha256']
            for cc in task['coefficients']:
                vv=tuple(map(F,cc));coeff_set.add(vv);expected,dd=make(r8,vv)
                assert expected==norm(task['input']) and dd==F(task['delta'])<=F(1,32)
                if index==0:bycoeff[vv]=rows[task['id']]
        if index==0:
            assert coeff_set=={(F(x,4),F(y,4),F(z,4)) for x,y,z in product(range(9),range(9),range(-4,5))}
        else:
            selection=m['selection'];center=tuple(map(F,selection['center']));axes=selection['axes'];scores=[]
            assert center==coarse_win
            win=bycoeff[center]
            for k in range(3):
                diffs=[]
                for sign in (-1,1):
                    vv=list(center);vv[k]+=sign*F(1,4)
                    if tuple(vv) in bycoeff:diffs.append(abs(F(win['inner_ratio'])-F(bycoeff[tuple(vv)]['inner_ratio'])))
                scores.append(max(diffs,default=F(0)))
            assert scores==list(map(F,selection['adjacent_coarse_inner_absolute_differences']))
            assert axes==sorted(range(3),key=lambda k:(scores[k],-k),reverse=True)[:2]==[0,1]
            expected_grid=set()
            for x,y in product(range(-16,17),repeat=2):
                vv=list(center);vv[axes[0]]+=F(x,64);vv[axes[1]]+=F(y,64);expected_grid.add(tuple(vv))
            assert coeff_set==expected_grid
        new=sum(not z.get('reused',False) for z in results);reuse=len(results)-new;hashes={z['input_sha256'] for z in results}
        assert new==row['new_calls'] and reuse==row['reused'] and len(hashes-all_hashes)==row['new_distinct_inputs']
        for z in results:
            assert 'error' not in z and z['input_sha256']==tasks[z['id']]['input_sha256']
            dd=F(tasks[z['id']]['delta']);bound=rho*(1+(F(r8['Px'])+2*dd)/10)*(1+(F(r8['Py'])+2*dd)/10)
            assert dd==F(z['delta']) and bound==F(z['tail_bound'])<1
            assert F(z['inner_ratio'])<=F(z['outer_ratio'])
            if not z.get('reused',False):actual[z['input_sha256']]+=1
        win=max(results,key=lambda z:F(z['inner_ratio']));assert win['id']==row['winner']['id']
        if index==0:coarse_win=tuple(map(F,tasks[win['id']]['coefficients'][0]));assert coarse_win==(F(1),F(1,2),F(0))
        else:assert tuple(map(F,tasks[win['id']]['coefficients'][0]))==tuple(coeff) and F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        all_hashes|=hashes;wall+=Decimal(str(row['wall_seconds']));cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        batches.append({'batch':index,'rows':len(results),'actual':new,'reused':reuse,'tasks':tasks,'results':results})
    for batch in batches:
        for z in batch['results']:
            if not z.get('reused',False):continue
            tasks,rows=filemap[z['reused_results_sha256']];prior=rows[z['reused_from_id']]
            assert z['input_sha256']==prior['input_sha256'] and F(z['delta'])==F(prior['delta'])
            assert F(z['inner_ratio'])==F(prior['inner_ratio']) and F(z['outer_ratio'])==F(prior['outer_ratio'])
            assert norm(batch['tasks'][z['id']]['input'])==norm(tasks[prior['id']]['input'])
    total=sum(z['rows'] for z in batches);reuse=sum(z['reused'] for z in batches)
    assert (total,len(all_hashes),sum(actual.values()),len(actual),reuse)==(1818,1809,1809,1809,9)
    assert f'{wall:.6f}'=='764.294603' and f'{cpu:.6f}'=='7600.528819'
    frac=lambda x:[F(int(u),int(v)) for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',x)]
    proof=t/'proofs/round014-lower.md';report=t/'rounds/014.md';assert frac(proof.read_text())==[tail,lo,up,gain]
    files=[ip,proof,report,res/'primary_final.json',res/'independent_final.json',res/'comparison.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','input_sha256':sha(ip),'R8_base_sha256':sha(r8p),'R12_baseline_sha256':sha(bp),
         'coefficients':list(map(str,coeff)),'positive_masses':204,'mass':str(mass),'changes_vs_R12':changes,
         'tail_radius':10,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R12_full_upper':str(gain),
         'mode_C_convention':'copy_y=0 group gets -C Py/100000 for copy_x=0 and +C Py/100000 for copy_x=1',
         'ledger_counts':{'rows':total,'distinct_inputs':len(all_hashes),'actual':sum(actual.values()),'reused':reuse,'duplicate_actual':sum(actual.values())-len(actual)},
         'fine_axes_and_scores_independently_verified':True,'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),
         'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':i['total_seconds'],
         'geometry_or_search_rerun':False,'files_sha256':{str(p):sha(p) for p in files},'audit_script_sha256':sha(Path(__file__)),
         'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
