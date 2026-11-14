"""Independent read-only reconstruction of every R17 task, LP alias and source difference."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from collections import Counter,defaultdict
from itertools import product
from datetime import datetime,timezone
import argparse,hashlib,json,math,re,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def norm(x):return {'Px':str(F(x['Px'])),'Py':str(F(x['Py'])),'lam':str(F(x['lam'])),'pts':[[str(F(v)) for v in q]for q in x['pts']]}
def weights(x):return list(map(lambda p:F(p[2]),x['pts']))
def mass_input(target,w):
    out=norm(target)
    for p,v in zip(out['pts'],w):p[2]=str(v)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round017'
    ip=t/'computation/inputs/r017_final_0015.json';b14=t/'computation/inputs/r014_final_0935.json';b15=t/'computation/inputs/r015_final_0680.json';b8=t/'computation/inputs/r008_final_0611.json'
    final,target,previous,base=map(read,(ip,b14,b15,b8));orig=weights(target);bw=weights(base);n=len(bw);assert n==51 and len(orig)==len(final['pts'])==204
    prior=raw.parent/'round015/supercell-mass';pm=read(prior/'batch0/manifest.json');pr=read(prior/'batch0/results.json');pt={x['id']:x for x in pm['tasks']};prows={x['id']:x for x in pr}
    ledger=read(r/'search_ledger.json');c=read(r/'integration_checks.json');p=read(r/'primary_final.json');ind=read(r/'independent_final.json')
    assert sha(ip)==p['input_sha256']==ind['input_sha256']==c['input_sha256']=='4ae1dd5c50d37fb6ffa05c75ab3c8e93e59e2b9beee0a13b09350442ff836595'
    assert sha(b14)==c['baseline_r14_sha256'] and sha(b15)==c['baseline_r15_sha256']
    assert all(F(u[k])==F(v[k]) for u,v in zip(final['pts'],target['pts']) for k in (0,1)) and all(F(final[k])==F(target[k])for k in ('Px','Py','lam'))
    px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));assert lam==1 and all(w>0 for w in weights(final))and mass==sum(orig)==4*sum(bw)
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts'])
    changed=[j for j,(x,y)in enumerate(zip(weights(previous),weights(final)))if x!=y];assert changed==c['mass_change_indices_from_r15']
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));gain=lo-F(read(t/'computation/results/round015/primary_final.json')['outer_ratio'])
    assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<1 and q==c['q']==3
    assert lo==F(ind['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(ind['weak_ratio_upper']['fraction'])==F(c['upper'])and gain==F(c['gain_over_r15_upper'])>0
    assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==mass*lo and F(p['outer_area'])==F(ind['area_upper']['fraction'])==mass*up
    assert p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==10
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'
    assert ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
    assert sha(raw/'mass-amplitudes/primary_old_R10.json')==sha(r/'primary_final.json') and sha(raw/'mass-amplitudes/candidate_0015.json')==sha(ip)
    # Reconstruct exact top eight ranked vectors directly from the prior complete coarse table.
    rank=sorted((z for z in pr if F(z['inner_ratio'])>F(prows['0000']['outer_ratio'])),key=lambda z:F(z['inner_ratio']),reverse=True);chosen=[];seen=set()
    for z in rank:
        task=pt[z['id']];v=tuple(x-y for x,y in zip(weights(task['input']),orig))
        if v in seen:continue
        seen.add(v);chosen.append({'source_id':z['id'],'source_input_sha256':z['input_sha256'],'aliases':task['aliases'],'vector':list(map(str,v)),'source_inner_ratio':z['inner_ratio']})
        if len(chosen)==8:break
    assert len(chosen)==8 and chosen==ledger[0]['selected_vectors'];vectors=[list(map(F,x['vector'])) for x in chosen]
    assert all(sum(v)==0 for v in vectors)
    P=[F(base['Px']),F(base['Py'])];delta=F(0)
    for row,point in enumerate(target['pts']):
        quad,j=divmod(row,n);cx,cy=divmod(quad,2)
        for axis,copy in ((0,cx),(1,cy)):
            d=(F(point[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(d,2*P[axis]-d))
    # Reconstruct all 128 LP columns and both finite-difference source pairs.
    lp=read(r/'mass_lp_input.json');num=read(r/'mass_lp_numeric.json');mo=json.loads((r/'mass_lp_output.json').read_text(),parse_float=Decimal)
    assert sha(Path(lp['base_input_path']))==lp['base_input_sha256']==sha(b8)and sha(Path(lp['target_input_path']))==lp['target_input_sha256']==sha(b14)
    assert lp['source_manifest_sha256']==sha(prior/'batch0/manifest.json')and lp['source_results_sha256']==sha(prior/'batch0/results.json')
    assert lp['export_script_sha256']==sha(t/'computation/export_mass_lp.py')and num['exact_input_sha256']==mo['exact_input_sha256']==sha(r/'mass_lp_input.json')
    assert list(map(F,lp['base_label_masses']))==bw*4 and list(map(F,lp['target_masses']))==orig
    low=list(map(F,lp['lower_mass_bound']));high=list(map(F,lp['upper_mass_bound']));assert low==[F(0)]*204 and high==[F(1001,1000)*w for w in bw]*4
    groups={}
    for task in pm['tasks']:
        for alias in task['aliases']:
            if alias.get('kind')=='baseline':continue
            key=(alias['quadrant'],tuple(alias['members']),alias['anchor']);g=groups.setdefault(key,{'aliases':[],'samples':{}});g['aliases'].append(alias)
            sk=(alias['eps'],alias['sign']);assert sk not in g['samples']or g['samples'][sk]==task['id'];g['samples'][sk]=task['id']
    assert len(groups)==len(lp['columns'])==128;cols=[];grads={s:[]for s in ('1/1000','1/10000')}
    for col,(key,g)in zip(lp['columns'],sorted(groups.items())):
        quad,members,anchor=key;assert (col['quadrant'],tuple(col['members']),col['anchor'])==key and col['aliases']==g['aliases']
        ids=[quad*n+j for j in members];aid=quad*n+anchor;G=sum(orig[j]for j in ids);cap=min(G,orig[aid]);v=[F(0)]*204
        for j in ids:v[j]=cap*orig[j]/G
        v[aid]-=cap;assert sum(v)==0 and list(map(F,col['unit_mass_increment']))==v;cols.append(v)
        for s in grads:
            e=col['estimates'][s];plus,minus=g['samples'][s,1],g['samples'][s,-1];prr,mrr=prows[plus],prows[minus]
            intervals=[[F(z[k])for k in ('inner_ratio','outer_ratio')]for z in (prr,mrr)]
            assert e['plus_id']==plus and e['minus_id']==minus and [list(map(F,e[k]))for k in ('plus_interval','minus_interval')]==intervals
            grad=(sum(intervals[0])-sum(intervals[1]))/(4*F(s));assert F(e['gradient_fraction_of_interval_midpoints'])==grad and e['gradient_float']==float(grad);grads[s].append(grad)
            # Also independently check each finite-difference input is +/- epsilon times this unit direction.
            for tid,sign in ((plus,1),(minus,-1)):
                assert weights(pt[tid]['input'])==[w+sign*F(s)*vv for w,vv in zip(orig,v)]
    assert num['W']==[[float(col[j])for col in cols]for j in range(204)]and num['G']==[[float(x)for x in grads[s]]for s in num['gradient_scales']]
    assert num['mass']==list(map(float,orig))and num['lower']==list(map(float,low))and num['upper']==list(map(float,high))
    assert num['budgets']==lp['budgets']==[1,2,4,8,16,32,64,132]and num['trust_exact']==lp['trust_scales']==['1/10000','3/10000','1/1000']and num['trust']==list(map(lambda s:float(F(s)),lp['trust_scales']))
    proposals=mo['proposals'];assert len(proposals)==mo['lp_count']==48 and all(z['exitflag']==1 for z in proposals)
    jobs={(z['gradient_index'],z['trust_index'],z['budget'])for z in proposals};assert jobs==set(product((1,2),(1,2,3),lp['budgets']))
    sparse=[[(j,x)for j,x in enumerate(v)if x]for v in cols];ray=[F(1,2),F(3,4),F(1),F(5,4),F(3,2),F(2)];expected_lp=defaultdict(list);alphas=Counter()
    for pi,z in enumerate(proposals):
        assert z['gradient_scale']==num['gradient_scales'][z['gradient_index']-1] and z['trust_exact']==lp['trust_scales'][z['trust_index']-1]and len(z['u'])==128
        assert all(x.is_finite() if isinstance(x,Decimal)else math.isfinite(x)for x in z['u'])
        # Independent round-to-even implementation uses exact Fraction.__round__ rather than Decimal rounding API.
        u=[F(round(F(x)*10**9),10**9)*F(z['trust_exact'])for x in z['u']]
        for scale in ray:
            dv=[F(0)]*204
            for coeff,col in zip(u,sparse):
                for j,x in col:dv[j]+=coeff*scale*x
            assert sum(dv)==0;alpha=F(1)
            for w,v,l,h in zip(orig,dv,low,high):
                if v>0:alpha=min(alpha,(h-w)/v)
                elif v<0:alpha=min(alpha,(w-l)/(-v))
            assert alpha>=0
            if any(w+alpha*v<=l for w,v,l in zip(orig,dv,low)):alpha*=F(999999999999,10**12)
            w=[x+alpha*v for x,v in zip(orig,dv)];assert sum(w)==sum(orig)and all(l<x<=h for x,l,h in zip(w,low,high))
            alias={'proposal_index':pi,'gradient_index':z['gradient_index'],'gradient_scale':z['gradient_scale'],'budget':z['budget'],'trust':z['trust_exact'],'ray_scale':str(scale),'exitflag':1,'feasibility_factor':str(alpha),'quantization':'nearest 1e-9 ties-to-even then exact trust'}
            expected_lp[digest(mass_input(target,w))].append(alias);alphas[str(alpha)]+=1
    assert sum(map(len,expected_lp.values()))==288 and len(expected_lp)==153
    # Inspect every saved task, result, completed row and every reuse source.
    filemap={}
    for d in (prior/'batch0',prior/'batch1'):
        mm=read(d/'manifest.json');rr=read(d/'results.json');filemap[sha(d/'results.json')]=({x['id']:x for x in mm['tasks']},{x['id']:x for x in rr})
    allhash=set();actual=Counter();stages=[];wall=Decimal(0);cpu=Decimal(0);reuse_rows=[]
    scans=['scan_binary_mass.py','scan_mass_amplitudes.py','scan_mass_lp.py'];modes=['binary-mass','mass-amplitudes','mass-lp/evaluation']
    for ix,mode in enumerate(modes):
        d=raw/mode;m=read(d/'manifest.json');rr=read(d/'results.json');done=[json.loads(x)for x in (d/'completed.jsonl').read_text().splitlines()];su=read(d/'summary.json');lr=ledger[ix]
        assert sha(d/'manifest.json')==lr['manifest_sha256']and sha(d/'results.json')==lr['results_sha256']and digest(m['tasks'])==m['task_list_sha256']==lr['task_list_sha256']
        assert m['scan_sha256']==lr['scan_sha256']==sha(t/'computation'/scans[ix])and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')
        assert m['bits']==40 and m['base_sha256']==sha(b8)and m['target_sha256']==sha(b14)and m['rejected']==[]and su['errors']==su['rejected']==0
        assert sorted(done,key=lambda x:x['id'])==sorted(rr,key=lambda x:x['id'])and len(rr)==len(m['tasks'])==su['completed']==lr['points']
        tasks={x['id']:x for x in m['tasks']};rows={x['id']:x for x in rr};assert len(tasks)==len(rows)==len(rr);filemap[sha(d/'results.json')]=(tasks,rows);combos=set()
        if ix<2:assert m['selected_vectors']==chosen
        if ix==0:assert m['prior_manifest_sha256']==sha(prior/'batch0/manifest.json')and m['prior_results_sha256']==sha(prior/'batch0/results.json')
        if ix==1:assert m['binary_manifest_sha256']==sha(raw/'binary-mass/manifest.json')
        if ix==2:assert m['lp_input_sha256']==sha(r/'mass_lp_input.json')and m['lp_output_sha256']==sha(r/'mass_lp_output.json')
        for tid,task in tasks.items():
            assert digest(task['input'])==task['input_sha256'];ww=weights(task['input']);assert sum(ww)==sum(orig)and all(w>0 for w in ww)
            assert norm(task['input'])==mass_input(target,ww)
            if ix<2:
                for alias in task['aliases']:
                    coef=tuple(alias['flags'if ix==0 else'coefficients']);combos.add(coef);rebuild=[w+sum(F(k)*v[j]for k,v in zip(coef,vectors))for j,w in enumerate(orig)];assert ww==rebuild
            else:assert task['aliases']==expected_lp[task['input_sha256']]
            row=rows[tid];assert 'error'not in row and row['input_sha256']==task['input_sha256'];k=max(w/bw[j%n]for j,w in enumerate(ww));R=F(row['tail_radius']);bound=k*sum(bw)/(P[0]*P[1])*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
            assert R==5 and k==F(row['kappa'])and delta==F(row['delta'])and bound==F(row['tail_bound'])<1 and F(row['inner_ratio'])<=F(row['outer_ratio'])
            if row.get('reused',False):reuse_rows.append((tasks,row))
            else:actual[row['input_sha256']]+=1
        if ix<2:assert combos==set(product((0,1)if ix==0 else(1,2),repeat=8))
        else:assert set(x['input_sha256']for x in tasks.values())==set(expected_lp)
        win=max(rr,key=lambda x:F(x['inner_ratio']));assert win['id']==lr['winner']['id'];hs={x['input_sha256']for x in rr};new=sum(not x.get('reused',False)for x in rr);reuse=len(rr)-new
        assert(new,reuse,len(hs-allhash))==(lr['new_calls'],lr['reused'],lr['new_distinct_inputs'])
        if ix==1:assert win['id']=='0015'and norm(tasks['0015']['input'])==norm(final)and F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
        if ix==2:
            lpbest=F(win['inner_ratio']);lpmaxouter=max(F(z['outer_ratio'])for z in rr);r14lo=F(read(t/'computation/results/round014/primary_final.json')['inner_ratio']);assert r14lo-lpbest==F(c['matlab_best_below_r14_lower'])>0 and r14lo-lpmaxouter==F(c['r14_lower_minus_matlab_max_outer'])>0
            assert c['matlab_max_outer_id']=='0000' and F(c['matlab_max_outer'])==lpmaxouter
        allhash|=hs;wall+=Decimal(str(lr['wall_seconds']));cpu+=Decimal(str(lr['children_user_cpu_seconds']))+Decimal(str(lr['children_system_cpu_seconds']));stages.append({'stage':mode,'rows':len(rr),'actual':new,'reuse':reuse,'aliases':sum(len(x['aliases'])for x in tasks.values())})
    for tasks,z in reuse_rows:
        oldtasks,oldrows=filemap[z['reused_results_sha256']];old=oldrows[z['reused_from_id']];assert z['input_sha256']==old['input_sha256']and norm(tasks[z['id']]['input'])==norm(oldtasks[old['id']]['input'])
        for k in ('kappa','delta','tail_bound','tail_radius','inner_ratio','outer_ratio'):assert F(z[k])==F(old[k])
    assert(sum(x['rows']for x in stages),len(allhash),sum(actual.values()),len(actual),len(reuse_rows))==(665,664,654,654,11)
    assert f'{wall:.6f}'=='284.253283'and f'{cpu:.6f}'=='2747.092103'
    coeff=[1,1,1,1,2,2,2,2];assert c['coefficients']==coeff and weights(final)==[w+sum(F(k)*v[j]for k,v in zip(coeff,vectors))for j,w in enumerate(orig)]
    for qd in range(4):
        for members,amp in (([32,33],1),([28,29],2)):
            for j in members:assert weights(final)[qd*n+j]==orig[qd*n+j]*(1+F(amp,10000))
        assert orig[qd*n+1]-weights(final)[qd*n+1]==sum(orig[qd*n+j]for j in(32,33))/10000+sum(orig[qd*n+j]for j in(28,29))/5000
    proof=t/'proofs/round017-lower.md';report=t/'rounds/017.md';fr=[F(int(u),int(v))for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]
    reportfr=[F(int(u),int(v))for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',report.read_text())];assert reportfr==[r14lo-lpmaxouter]
    for short,full in [('input','mass_lp_input.json'),('numeric','mass_lp_numeric.json'),('output','mass_lp_output.json')]:assert sha(raw/'mass-lp'/f'{short}.json')==sha(r/full)
    files=[ip,proof,report,*[t/'computation'/x for x in scans+['export_mass_lp.py','build_mass_lp_numeric.py','propose_mass_joint.m']],*sorted(r.glob('*.json'))]
    out={'status':'passed','input_sha256':sha(ip),'positive_masses':204,'mass':str(mass),'periods':[str(px),str(py)],'positions_unchanged_from_R14':True,'changed_mass_indices_from_R15':changed,'winner_coefficients':coeff,'selected_source_ids':[x['source_id']for x in chosen],
         'tail_radius':10,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R15_full_upper':str(gain),'stages':stages,'ledger_counts':{'rows':665,'unique_inputs':664,'actual':654,'unique_actual':654,'reuse':11},
         'lp_columns_reconstructed':128,'finite_difference_source_pairs_reconstructed':256,'formal_matlab_proposals':48,'matlab_exitflags':{'1':48},'matlab_version':mo['version'],'matlab_recorded_wall_seconds':str(mo['wall_seconds']),'lp_aliases_reconstructed':288,'lp_unique_candidates_reconstructed':153,'feasibility_factors':dict(alphas),'LP_best_inner_below_R14_inner':str(r14lo-lpbest),'LP_maximum_outer':str(lpmaxouter),'R14_inner_minus_LP_maximum_outer':str(r14lo-lpmaxouter),
         'heuristic_scope':'Finite differences are midpoint estimates, not derivatives; exact ray repair certifies mass constraints, not numerical LP optimality or trust budget.','search_wall_sum':str(wall),'child_cpu_sum':str(cpu),'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'geometry_or_search_rerun':False,'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'python':sys.version,'elapsed_seconds':time.perf_counter()-start,'verified_utc':datetime.now(timezone.utc).isoformat(),'command':[sys.executable,*sys.argv]}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in('files_sha256','feasibility_factors')},indent=2))
if __name__=='__main__':main()
