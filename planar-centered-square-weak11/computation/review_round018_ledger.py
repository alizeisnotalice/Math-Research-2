"""Independent exact response-selection and complete joint-grid ledger audit for R18."""
from review_round017_ledger import *

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round018';root=raw.parent.parent
    ip=t/'computation/inputs/r018_final_0728.json';prevp=t/'computation/inputs/r017_final_0015.json';b14=t/'computation/inputs/r014_final_0935.json';b8=t/'computation/inputs/r008_final_0611.json'
    final,previous,target,base=map(read,(ip,prevp,b14,b8));orig=weights(target);bw=weights(base);n=len(bw);assert n==51 and len(final['pts'])==len(orig)==204
    p,ind,c=map(read,(r/'primary_final.json',r/'independent_final.json',r/'integration_checks.json'));ledger=read(r/'search_ledger.json')
    assert sha(ip)==p['input_sha256']==ind['input_sha256']==c['input_sha256']=='41e4cc6ff3d0d0dff12c04dc280ad64fd203d0b4278be3ecfb0273da53ca7726'and sha(prevp)==c['baseline_r17_sha256']
    assert all(F(u[k])==F(v[k])for u,v in zip(final['pts'],previous['pts'])for k in (0,1))and all(F(final[k])==F(previous[k])for k in('Px','Py','lam'))
    px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));assert lam==1 and all(w>0 for w in weights(final))and mass==sum(orig)==sum(weights(previous))==4*sum(bw)
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts']);changed=[j for j,(x,y)in enumerate(zip(weights(previous),weights(final)))if x!=y];assert changed==c['mass_change_indices_from_r17']
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));gain=lo-F(read(t/'computation/results/round017/primary_final.json')['outer_ratio'])
    assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    assert lo==F(ind['weak_ratio_lower']['fraction'])==F(c['lower'])and up==F(ind['weak_ratio_upper']['fraction'])==F(c['upper'])and gain==F(c['gain_over_r17_upper'])>0
    assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(ind['area_upper']['fraction'])==up*mass and F(p['tail_radius'])==F(ind['radius_cap'])==10 and p['bits']==ind['bits']==40
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
    # Raw certificate byte equality and unchanged exact source-vector evidence.
    assert sha(raw/'primary_old_R10.json')==sha(r/'primary_final.json')and sha(raw/'batch1/candidate_0728.json')==sha(ip)
    prior=root/'round017';bm=read(prior/'binary-mass/manifest.json');am=read(prior/'mass-amplitudes/manifest.json');vectors=[list(map(F,x['vector']))for x in bm['selected_vectors']]
    assert am['selected_vectors']==bm['selected_vectors']and all(sum(v)==0 for v in vectors)and bm['target_sha256']==am['target_sha256']==sha(b14)
    audit17=read(prior/'independent-ledger-review.json');assert audit17['status']=='passed'and audit17['input_sha256']==sha(prevp)
    assert audit17['files_sha256'][str(t/'computation/results/round017/search_ledger.json')]==sha(t/'computation/results/round017/search_ledger.json')
    amp_summary=read(prior/'mass-amplitudes/summary.json');center=tuple(map(F,amp_summary['top'][0]['aliases'][0]['coefficients']));assert center==(1,1,1,1,2,2,2,2)
    def rebuild(coef):return mass_input(target,[w+sum(k*v[j]for k,v in zip(coef,vectors))for j,w in enumerate(orig)])
    assert norm(previous)==rebuild(center)
    P=[F(base['Px']),F(base['Py'])];delta=F(0)
    for row,point in enumerate(target['pts']):
        quad,j=divmod(row,n);cx,cy=divmod(quad,2)
        for axis,copy in((0,cx),(1,cy)):
            z=(F(point[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(z,2*P[axis]-z))
    filemap={};cachehash=set()
    for d in [root/'round015/supercell-mass/batch0',root/'round015/supercell-mass/batch1',prior/'binary-mass',prior/'mass-amplitudes']:
        m=read(d/'manifest.json');rr=read(d/'results.json');filemap[sha(d/'results.json')]=({x['id']:x for x in m['tasks']},{x['id']:x for x in rr})
        if d.parent==prior:cachehash|={x['input_sha256']for x in rr}
    steps=[F(-1),F(-1,2),F(1,2),F(1)];expected0={center}
    for j in range(8):
        for step in steps:
            v=list(center);v[j]+=step;expected0.add(tuple(v))
    stages=[];allhash=set();actual=Counter();reuse_records=[];wall=Decimal(0);cpu=Decimal(0);directions=None;batchhash=[]
    for ix,lr in enumerate(ledger):
        d=raw/f'batch{ix}';m=read(d/'manifest.json');rr=read(d/'results.json');done=[json.loads(x)for x in(d/'completed.jsonl').read_text().splitlines()];su=read(d/'summary.json')
        assert sha(d/'manifest.json')==lr['manifest_sha256']and sha(d/'results.json')==lr['results_sha256']and digest(m['tasks'])==m['task_list_sha256']==lr['task_list_sha256']
        assert m['scan_sha256']==lr['scan_sha256']==sha(t/'computation/scan_mass_response.py')and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')and m['base_sha256']==sha(b8)and m['target_sha256']==sha(b14)
        assert m['selection']==lr['selection']and m['bits']==40 and m['rejected']==[]and su['errors']==su['rejected']==0 and F(m['delta'])==delta
        assert len(rr)==len(m['tasks'])==len(done)==su['completed']==lr['points']and sorted(done,key=lambda x:x['id'])==sorted(rr,key=lambda x:x['id'])
        tasks={x['id']:x for x in m['tasks']};rows={x['id']:x for x in rr};assert len(tasks)==len(rows)==len(rr);filemap[sha(d/'results.json')]=(tasks,rows);bycoef={}
        for tid,task in tasks.items():
            assert digest(task['input'])==task['input_sha256']
            for co in task['coefficients']:
                coef=tuple(map(F,co));assert coef not in bycoef;bycoef[coef]=rows[tid];assert norm(task['input'])==rebuild(coef)
            ww=weights(task['input']);assert sum(ww)==sum(orig)and min(ww)>0
            z=rows[tid];assert 'error'not in z and z['input_sha256']==task['input_sha256'];k=max(w/bw[j%n]for j,w in enumerate(ww));R=F(z['tail_radius']);bound=k*sum(bw)/(P[0]*P[1])*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
            assert R==5 and k==F(z['kappa'])and F(z['delta'])==delta and bound==F(z['tail_bound'])<1 and F(z['inner_ratio'])<=F(z['outer_ratio'])
            if z.get('reused',False):reuse_records.append((tasks,z))
            else:actual[z['input_sha256']]+=1
        if ix==0:
            assert set(bycoef)==expected0 and m['selection']['source_binary_sha256']==sha(prior/'binary-mass/manifest.json')and m['selection']['source_amplitudes_sha256']==sha(prior/'mass-amplitudes/manifest.json')
            baseline=bycoef[center];directions=[]
            for j in range(8):
                eligible=[]
                for step in steps:
                    cc=list(center);cc[j]+=step;z=bycoef[tuple(cc)]
                    if F(z['inner_ratio'])>F(baseline['outer_ratio']):eligible.append((F(z['inner_ratio']),step,z['id']))
                if eligible:
                    score,step,tid=max(eligible);directions.append((score,j,step,tid))
            directions=sorted(directions,reverse=True)[:6];assert len(directions)==6
        else:
            desc=[{'index':j,'step':str(step),'source_id':tid,'source_lower':str(score)}for score,j,step,tid in directions];assert m['selection']['directions']==desc
            expected=set()
            for signs in product((-1,0,1),repeat=6):
                cc=list(center)
                for sign,(_,j,step,_)in zip(signs,directions):cc[j]+=sign*step
                expected.add(tuple(cc))
            assert set(bycoef)==expected and len(expected)==729
            win=max(rr,key=lambda z:F(z['inner_ratio']));assert win['id']==lr['winner']['id']=='0728'and norm(tasks[win['id']]['input'])==norm(final)and F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        hs={z['input_sha256']for z in rr};new=sum(not z.get('reused',False)for z in rr);reuse=len(rr)-new;assert(new,reuse,len(hs-allhash))==(lr['new_calls'],lr['reused'],lr['new_distinct_inputs'])
        batchhash.append(hs);allhash|=hs;wall+=Decimal(str(lr['wall_seconds']));cpu+=Decimal(str(lr['children_user_cpu_seconds']))+Decimal(str(lr['children_system_cpu_seconds']));stages.append({'batch':ix,'rows':len(rr),'actual':new,'reused':reuse,'all_aliases':len(bycoef)})
    checked_sources=set()
    def trace(tasks,z):
        if not z.get('reused',False):return
        key=z['reused_results_sha256'],z['reused_from_id'];oldtasks,oldrows=filemap[key[0]];old=oldrows[key[1]];assert z['input_sha256']==old['input_sha256']and norm(tasks[z['id']]['input'])==norm(oldtasks[old['id']]['input'])
        for k in('kappa','delta','tail_radius','tail_bound','inner_ratio','outer_ratio'):assert F(z[k])==F(old[k])
        if key not in checked_sources:checked_sources.add(key);trace(oldtasks,old)
    for tasks,z in reuse_records:trace(tasks,z)
    assert(sum(x['rows']for x in stages),len(allhash),sum(actual.values()),len(actual),len(reuse_records),len(batchhash[0]&batchhash[1]),len(allhash-set(actual)))==(762,749,729,729,33,13,20)
    assert f'{wall:.6f}'=='301.544018'and f'{cpu:.6f}'=='3043.112718'
    coeff=list(map(F,c['coefficients']));assert coeff==list(map(F,['1/2','1/2','1','1','3','3','3','3']))and norm(final)==rebuild(coeff)
    proof=t/'proofs/round018-lower.md';report=t/'rounds/018.md';fr=[F(int(u),int(v))for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]
    files=[ip,proof,report,t/'computation/scan_mass_response.py',*sorted(r.glob('*.json'))]
    out={'status':'passed','input_sha256':sha(ip),'positive_masses':204,'mass':str(mass),'periods':[str(px),str(py)],'coordinates_and_periods_unchanged':True,'changed_mass_indices_from_R17':changed,'winner_coefficients':list(map(str,coeff)),'selected_six_directions':[{'axis':j,'step':str(step),'source_id':tid,'source_lower':str(score)}for score,j,step,tid in directions],
         'tail_radius':10,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R17_full_upper':str(gain),'stages':stages,'ledger_counts':{'rows':762,'unique_inputs':749,'actual':729,'unique_actual':729,'reuse':33,'shared_between_batches':13,'distinct_from_earlier_batches':20},'wall_sum':str(wall),'child_cpu_sum':str(cpu),'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'geometry_or_search_rerun':False,'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'python':sys.version,'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
