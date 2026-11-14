"""Read-only exact R12 construction, saved old-tail certificate and ledger audit."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from collections import Counter
from datetime import datetime,timezone
import argparse,hashlib,json,math,re,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def norm(x):return {'Px':str(F(x['Px'])),'Py':str(F(x['Py'])),'lam':str(F(x['lam'])),
                    'pts':[[str(F(v)) for v in q] for q in x['pts']]}
def make(base,alias):
    px,py=F(base['Px']),F(base['Py']);pts=[];delta=F(0);group=set(range(27,50,2))
    axis=alias.get('axis');members=alias.get('members',[]);shift=F(alias.get('shift','0'))
    if 'coordinate' in alias:
        assert members==[i for i,p in enumerate(base['pts']) if F(p[axis])==F(alias['coordinate'])]
    for a in (0,1):
        for b in (0,1):
            for i,p in enumerate(base['pts']):
                x,y,w=map(F,p);hx=F(0);hy=3*py/400000 if b==1 and i in group else F(0)
                if a==b==1 and i in members:
                    if axis==0:hx+=shift*px
                    else:hy+=shift*py
                delta=max(delta,abs(hx),abs(hy))
                pts.append([str((x+a*px+hx)%(2*px)),str((y+b*py+hy)%(2*py)),str(w)])
    return {'Px':str(2*px),'Py':str(2*py),'lam':str(F(base['lam'])),'pts':pts},delta

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();start=time.perf_counter();t=a.topic;res=t/'computation/results/round012'
    ip=t/'computation/inputs/r012_final_0017.json';bp=t/'computation/inputs/r010_final_0019.json'
    r8p=t/'computation/inputs/r008_final_0611.json';d,b,r8=map(read,(ip,bp,r8p))
    p,i,c=map(read,(res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json'))
    old=read(t/'computation/results/round010/primary_final.json');ledger=read(res/'search_ledger.json')
    assert sha(ip)==p['input_sha256']==i['input_sha256']==c['input_sha256'] and sha(bp)==c['baseline_sha256']
    assert sha(a.raw_search/'primary_old_R10.json')==sha(res/'primary_final.json')
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py') and i['script_sha256']==sha(t/'computation/independent_periodic_halo.py')
    assert F(p['tail_radius'])==F(i['radius_cap'])==10 and p['bits']==i['bits']==40
    expected,delta=make(r8,{'kind':'quadrant_shift','axis':1,'members':[50],'shift':'1/400000'})
    assert expected==norm(d)
    baseline,_=make(r8,{'kind':'baseline'})
    px0,py0=F(r8['Px']),F(r8['Py'])
    # Verify the reference is exactly two horizontal copies of frozen R10.
    expected_baseline=[]
    for copy in (0,1):
        for x,y,w in b['pts']:expected_baseline.append([str((F(x)+copy*F(b['Px']))%(2*F(b['Px']))),str(F(y)),str(F(w))])
    assert baseline['pts']==expected_baseline
    diffs=[(j,k,str(F(z)-F(y))) for j,(p0,p1) in enumerate(zip(expected_baseline,d['pts'])) for k,(y,z) in enumerate(zip(p0,p1)) if F(y)!=F(z)]
    assert diffs==[(203,1,str(py0/400000))]
    mass=sum(F(q[2]) for q in d['pts']);mass10=sum(F(q[2]) for q in b['pts']);mass8=sum(F(q[2]) for q in r8['pts'])
    px,py,lam=map(F,(d['Px'],d['Py'],d['lam']))
    assert len(d['pts'])==204 and mass==2*mass10==4*mass8
    assert (px,py,lam)==(2*px0,2*py0,F(1)) and all(F(p[2])>0 for p in d['pts'])
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1
    assert tail==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    lo,up=F(p['inner_ratio']),F(p['outer_ratio']);margin=lo-F(old['outer_ratio'])
    assert lo==F(i['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(i['weak_ratio_upper']['fraction'])==F(c['upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lo*mass
    assert F(p['outer_area'])==F(i['area_upper']['fraction'])==up*mass and margin==F(c['gain_over_baseline_upper'])>0
    assert i['counts']['all_span_admissible_boxes']==42486322
    assert i['counts']['unique_inner_rectangles']==i['counts']['unique_outer_rectangles']==23904
    assert F(c['epsilon'])==F(1,2**41)
    all_hashes=set();actual=Counter();filemap={};batches=[];total=reused=0;wall=Decimal(0);cpu=Decimal(0)
    rho=mass8/(px0*py0)
    for index,row in enumerate(ledger):
        raw=a.raw_search/f'batch{index}';m=read(raw/'manifest.json');results=read(raw/'results.json')
        completed=[json.loads(x) for x in (raw/'completed.jsonl').read_text().splitlines()]
        assert sha(raw/'manifest.json')==row['manifest_sha256'] and sha(raw/'results.json')==row['results_sha256']
        assert digest(m['tasks'])==row['task_list_sha256']==m['task_list_sha256']
        assert m['input_sha256']==sha(r8p) and m['copy_factors']==[2,2] and m['radius']==5
        assert m['transport_kernel_sha256']==row['transport_kernel_sha256']==sha(t/'computation/periodic_halo_transport.py')
        assert len(results)==len(completed)==len(m['tasks'])==row['points']
        assert sorted(completed,key=lambda z:z['id'])==sorted(results,key=lambda z:z['id'])
        assert len({z['id'] for z in results})==len(results)
        tasks={z['id']:z for z in m['tasks']};byid={z['id']:z for z in results};filemap[row['results_sha256']]=(tasks,byid)
        for task in tasks.values():
            assert digest(task['input'])==task['input_sha256']
            for alias in task['aliases']:
                expected,dd=make(r8,alias)
                assert expected==norm(task['input']) and dd==F(task['delta'])<=F(1,32)
        hashes={z['input_sha256'] for z in results};nr=sum(not z.get('reused',False) for z in results);rr=len(results)-nr
        assert nr==row['new_worker_evaluations'] and rr==row['reused'] and len(hashes-all_hashes)==row['new_distinct_inputs']
        for z in results:
            assert 'error' not in z and z['input_sha256']==tasks[z['id']]['input_sha256']
            dd=F(tasks[z['id']]['delta']);g=rho*(1+(px0+2*dd)/10)*(1+(py0+2*dd)/10)
            assert dd==F(z['delta']) and g==F(z['tail_bound'])<1
            assert F(z['inner_ratio'])<=F(z['outer_ratio'])
            if not z.get('reused',False):actual[z['input_sha256']]+=1
        win=max(results,key=lambda z:F(z['inner_ratio']));assert win['id']==row['winner']['id']
        assert win['input_sha256']==row['winner']['input_sha256']
        if index==1:assert F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        total+=len(results);reused+=rr;all_hashes|=hashes
        wall+=Decimal(str(row['wall_seconds']));cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        batches.append({'batch':index,'rows':len(results),'actual':nr,'reused':rr,'tasks':tasks,'results':results})
    for batch in batches:
        for z in batch['results']:
            if not z.get('reused',False):continue
            tasks,rows=filemap[z['reused_results_sha256']];oldrow=rows[z['reused_from_id']]
            assert z['input_sha256']==oldrow['input_sha256']
            assert norm(batch['tasks'][z['id']]['input'])==norm(tasks[oldrow['id']]['input'])
            assert F(z['delta'])==F(oldrow['delta']) and F(z['inner_ratio'])==F(oldrow['inner_ratio']) and F(z['outer_ratio'])==F(oldrow['outer_ratio'])
    assert (total,len(all_hashes),sum(actual.values()),len(actual),reused)==(166,163,163,163,3)
    assert f'{wall:.6f}'=='61.248662' and f'{cpu:.6f}'=='599.411406'
    frac=lambda x:[F(int(a),int(b)) for a,b in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',x)]
    proof=t/'proofs/round012-lower.md';report=t/'rounds/012.md';assert frac(proof.read_text())==[tail,lo,up,margin]
    r1=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/periodic-complete-halo.md')
    assert sha(r1)=='7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e'
    files=[ip,proof,report,res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','input_sha256':sha(ip),'baseline_R10_sha256':sha(bp),'positive_masses':204,'mass':str(mass),
         'extra_changed_label':[1,1,50],'extra_y_shift':str(py0/400000),'old_tail_radius':10,'tail_bound':str(tail),'q':q,
         'lower':str(lo),'upper':str(up),'lower_minus_R10_full_upper':str(margin),
         'primary_reused_file_sha256':sha(a.raw_search/'primary_old_R10.json'),
         'independent_recorded_seconds':i['total_seconds'],'independent_enumerated_boxes':42486322,
         'ledger_counts':{'rows':total,'distinct_inputs':len(all_hashes),'actual_search_evaluations':sum(actual.values()),'reuse':reused,'duplicate_actual_evaluations':sum(actual.values())-len(actual)},
         'search_uses_single_R5_transport_call_per_actual_input':True,
         'batches':[{'batch':z['batch'],'rows':z['rows'],'actual':z['actual'],'reused':z['reused']} for z in batches],
         'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),'geometry_or_search_rerun':False,
         'files_sha256':{str(p):sha(p) for p in files},'audit_script_sha256':sha(Path(__file__)),
         'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
