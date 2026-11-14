"""Independent exact audit of the complete R20 two-parameter grids and saved full certificates."""
from review_round017_ledger import *

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round020';root=raw.parent.parent
    ip=t/'computation/inputs/r020_final_0730.json';prevp=t/'computation/inputs/r018_final_0728.json';b14=t/'computation/inputs/r014_final_0935.json';b8=t/'computation/inputs/r008_final_0611.json'
    final,previous,target,base=map(read,(ip,prevp,b14,b8));orig=weights(target);bw=weights(base);n=len(bw);assert n==51 and len(final['pts'])==len(orig)==204
    p,ind,c=map(read,(r/'primary_final.json',r/'independent_final.json',r/'integration_checks.json'));ledger=read(r/'search_ledger.json')
    assert sha(ip)==p['input_sha256']==ind['input_sha256']==c['input_sha256']=='b8fa96aa7dad1c2fa0f2927e6e98049fc42be132897c950662fcab3fcbe9df5a'and sha(prevp)==c['baseline_r18_sha256']
    assert all(F(u[k])==F(v[k])for u,v in zip(final['pts'],previous['pts'])for k in(0,1))and all(F(final[k])==F(previous[k])for k in('Px','Py','lam'))
    px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));assert lam==1 and min(weights(final))>0 and mass==sum(orig)==sum(weights(previous))==4*sum(bw)
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts']);changed=[j for j,(x,y)in enumerate(zip(weights(previous),weights(final)))if x!=y];assert changed==c['mass_change_indices_from_r18']
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));gain=lo-F(read(t/'computation/results/round018/primary_final.json')['outer_ratio'])
    assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    assert lo==F(ind['weak_ratio_lower']['fraction'])==F(c['lower'])and up==F(ind['weak_ratio_upper']['fraction'])==F(c['upper'])and gain==F(c['gain_over_r18_upper'])>0
    assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(ind['area_upper']['fraction'])==up*mass and F(p['tail_radius'])==F(ind['radius_cap'])==10 and p['bits']==ind['bits']==40
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
    assert sha(raw/'primary_old_R10.json')==sha(r/'primary_final.json')and sha(raw/'batch1/candidate_0730.json')==sha(ip)
    # Source family was independently fully reconstructed in R17 and R18.
    prior=root/'round017';bm=read(prior/'binary-mass/manifest.json');vectors=[list(map(F,x['vector']))for x in bm['selected_vectors']];assert len(vectors)==8 and all(sum(v)==0 for v in vectors)
    audit17=read(prior/'independent-ledger-review.json');audit18=read(root/'round018/independent-ledger-review.json');assert audit17['status']==audit18['status']=='passed'and audit18['input_sha256']==sha(prevp)
    assert audit17['files_sha256'][str(t/'computation/results/round017/search_ledger.json')]==sha(t/'computation/results/round017/search_ledger.json')
    def convert(x,y):return (x,x,F(1),F(1),y,y,y,y)
    def rebuild(coef):return mass_input(target,[w+sum(k*v[j]for k,v in zip(coef,vectors))for j,w in enumerate(orig)])
    aa,bb=map(F,c['parameters']);assert (aa,bb)==(F(-3,20),F(33,4))and norm(final)==rebuild(convert(aa,bb))
    high=[F(1001,1000)*w for w in bw]*4;P=[F(base['Px']),F(base['Py'])];delta=F(0)
    for row,point in enumerate(target['pts']):
        quad,j=divmod(row,n);cx,cy=divmod(quad,2)
        for axis,copy in((0,cx),(1,cy)):
            z=(F(point[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(z,2*P[axis]-z))
    provenance=read(raw/'orchestrator_provenance.json');orchestrator=t/'computation/scan_mass_common.py';runner=t/'computation/scan_mass_response.py'
    assert provenance['orchestrator']=='scan_mass_common.py'and provenance['orchestrator_sha256']==sha(orchestrator)and provenance['shared_batch_worker']=='scan_mass_response.py'and provenance['shared_batch_worker_sha256']==sha(runner)
    if (r/'orchestrator_provenance.json').exists():assert sha(r/'orchestrator_provenance.json')==sha(raw/'orchestrator_provenance.json')
    filemap={}
    for d in [root/'round015/supercell-mass/batch0',root/'round015/supercell-mass/batch1',prior/'binary-mass',prior/'mass-amplitudes',root/'round018/mass-response/batch0',root/'round018/mass-response/batch1']:
        m=read(d/'manifest.json');rr=read(d/'results.json');filemap[sha(d/'results.json')]=({x['id']:x for x in m['tasks']},{x['id']:x for x in rr})
    # Exact Cartesian products, including the separately appended R18 reference point.
    expected0={convert(F(i,8),F(5*j,8))for i,j in product(range(-8,9),range(17))}|{convert(F(1,2),F(3))}
    expected1={convert(F(i,20),F(j,4))for i,j in product(range(-20,21),range(41))};assert len(expected0)==290 and len(expected1)==1681
    stages=[];allhash=set();actual=Counter();reuse_records=[];wall=Decimal(0);cpu=Decimal(0);batchhash=[];coarse_win=None
    for ix,lr in enumerate(ledger):
        d=raw/f'batch{ix}';m=read(d/'manifest.json');rr=read(d/'results.json');done=[json.loads(x)for x in(d/'completed.jsonl').read_text().splitlines()];su=read(d/'summary.json')
        assert sha(d/'manifest.json')==lr['manifest_sha256']and sha(d/'results.json')==lr['results_sha256']and digest(m['tasks'])==m['task_list_sha256']==lr['task_list_sha256']
        assert m['scan_sha256']==lr['run_kernel_sha256']==sha(runner)and lr['grid_source_sha256']==sha(orchestrator)and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')and m['base_sha256']==sha(b8)and m['target_sha256']==sha(b14)
        assert m['selection']==lr['selection']and m['bits']==40 and m['rejected']==m['selection']['excluded']==[]and m['selection']['kappa_cap']=='1001/1000'and su['errors']==su['rejected']==0 and F(m['delta'])==delta
        assert len(rr)==len(m['tasks'])==len(done)==su['completed']==lr['points']and sorted(done,key=lambda x:x['id'])==sorted(rr,key=lambda x:x['id'])
        tasks={x['id']:x for x in m['tasks']};rows={x['id']:x for x in rr};assert len(tasks)==len(rows)==len(rr);filemap[sha(d/'results.json')]=(tasks,rows);bycoef={}
        for tid,task in tasks.items():
            assert digest(task['input'])==task['input_sha256']
            for co in task['coefficients']:
                coef=tuple(map(F,co));assert coef not in bycoef;bycoef[coef]=rows[tid];assert norm(task['input'])==rebuild(coef)
            ww=weights(task['input']);assert sum(ww)==sum(orig)and min(ww)>0 and all(w<=h for w,h in zip(ww,high))
            z=rows[tid];assert 'error'not in z and z['input_sha256']==task['input_sha256'];k=max(w/bw[j%n]for j,w in enumerate(ww));R=F(z['tail_radius']);bound=k*sum(bw)/(P[0]*P[1])*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
            assert R==5 and k==F(z['kappa'])<=F(1001,1000)and F(z['delta'])==delta and bound==F(z['tail_bound'])<1 and F(z['inner_ratio'])<=F(z['outer_ratio'])
            if z.get('reused',False):reuse_records.append((tasks,z))
            else:actual[z['input_sha256']]+=1
        assert set(bycoef)==(expected0 if ix==0 else expected1)
        win=max(rr,key=lambda z:F(z['inner_ratio']));assert win['id']==lr['winner']['id'];wc=tuple(map(F,tasks[win['id']]['coefficients'][0]))
        if ix==0:
            assert m['selection']['source_binary_sha256']==sha(prior/'binary-mass/manifest.json');coarse_win=[str(wc[0]),str(wc[4])];assert coarse_win==['-1/8','65/8']
        else:
            assert m['selection']['coarse_winner']==coarse_win and win['id']=='0730'and norm(tasks[win['id']]['input'])==norm(final)and F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        hs={z['input_sha256']for z in rr};new=sum(not z.get('reused',False)for z in rr);reuse=len(rr)-new;assert(new,reuse,len(hs-allhash))==(lr['new_calls'],lr['reused'],lr['new_distinct_inputs'])
        batchhash.append(hs);allhash|=hs;wall+=Decimal(str(lr['wall_seconds']));cpu+=Decimal(str(lr['children_user_cpu_seconds']))+Decimal(str(lr['children_system_cpu_seconds']));stages.append({'batch':ix,'rows':len(rr),'actual':new,'reused':reuse,'all_aliases':len(bycoef)})
    checked_sources=set()
    def trace(tasks,z):
        if not z.get('reused',False):return
        key=z['reused_results_sha256'],z['reused_from_id'];oldtasks,oldrows=filemap[key[0]];old=oldrows[key[1]];assert z['input_sha256']==old['input_sha256']and norm(tasks[z['id']]['input'])==norm(oldtasks[old['id']]['input'])
        for k in('kappa','delta','tail_radius','tail_bound','inner_ratio','outer_ratio'):assert F(z[k])==F(old[k])
        if key not in checked_sources:checked_sources.add(key);trace(oldtasks,old)
    for tasks,z in reuse_records:trace(tasks,z)
    assert(sum(x['rows']for x in stages),len(allhash),sum(actual.values()),len(actual),len(reuse_records),len(batchhash[0]&batchhash[1]),len(allhash-set(actual)))==(1971,1889,1880,1880,91,82,9)
    assert f'{wall:.6f}'=='797.818018'and f'{cpu:.6f}'=='7749.773747'
    proof=t/'proofs/round020-lower.md';report=t/'rounds/020.md';fr=[F(int(u),int(v))for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]
    files=[ip,proof,report,orchestrator,runner,raw/'orchestrator_provenance.json',*sorted(r.glob('*.json'))]
    out={'status':'passed','input_sha256':sha(ip),'positive_masses':204,'mass':str(mass),'periods':[str(px),str(py)],'coordinates_and_periods_unchanged':True,'changed_mass_indices_from_R18':changed,'winner_parameters':[str(aa),str(bb)],'winner_coefficients':list(map(str,convert(aa,bb))),'tail_radius':10,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R18_full_upper':str(gain),'stages':stages,'ledger_counts':{'rows':1971,'unique_inputs':1889,'actual':1880,'unique_actual':1880,'reuse':91,'shared_between_batches':82,'distinct_from_earlier_batches':9},'grid_orchestrator_sha256':sha(orchestrator),'shared_run_worker_sha256':sha(runner),'provenance_archived':(r/'orchestrator_provenance.json').exists(),'wall_sum':str(wall),'child_cpu_sum':str(cpu),'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'geometry_or_search_rerun':False,'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'audit_helper_sha256':sha(Path(__file__).with_name('review_round017_ledger.py')),'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
