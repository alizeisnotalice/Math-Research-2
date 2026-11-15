"""Independent exact verification of R26 coordinate groups, reuse, and C2 lower bound."""
from review_round017_ledger import *
import gzip

def displacement(base,target):
    P=[F(base['Px']),F(base['Py'])];delta=F(0);n=len(base['pts'])
    for k,row in enumerate(target['pts']):
        cx,cy=divmod(k//n,2);j=k%n
        for axis,copy in ((0,cx),(1,cy)):
            d=(F(row[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(d,2*P[axis]-d))
    return delta

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round026';prevraw=raw.parent.parent/'round025/fixedmass-position'
    bp=t/'computation/inputs/r008_final_0611.json';oldp=t/'computation/inputs/r025_final_0110.json';ip=t/'computation/inputs/r026_final_0360.json';base,old,final=map(read,(bp,oldp,ip));m=read(raw/'manifest.json');rows=read(raw/'results.json');su=read(raw/'summary.json');done=[json.loads(s)for s in(raw/'completed.jsonl').read_text().splitlines()];diag=read(raw/'geometry_diagnostics.json');c=read(r/'integration_checks.json')
    assert sha(ip)=='eaa5e71e4b5060eb804208fd2fa853e5cff4e5a7856ac9b78f7210bf85494aa0'and sha(oldp)=='34a9569786741ddf2a61cc52fd85725bf8efdccf46306e34c1d8d340be49cd0c'
    n=len(base['pts']);P=[F(base['Px']),F(base['Py'])];px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));assert n==51 and len(old['pts'])==len(final['pts'])==204 and min(weights(final))>0 and weights(final)==weights(old)and mass==4*sum(weights(base))and(px,py)==tuple(2*v for v in P)and lam==1
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts'])
    changes={name:[j for j,(u,v)in enumerate(zip(old['pts'],final['pts']))if F(u[k])!=F(v[k])]for k,name in enumerate(('x','y','weight'))};assert changes=={'x':[79,80],'y':[],'weight':[]}
    for k in (79,80):assert F(final['pts'][k][0])-F(old['pts'][k][0])==P[0]/500000
    assert all(F(final[k])==F(old[k])for k in('Px','Py','lam'))
    # Reconstruct every original coordinate equivalence class, then all signed moves.
    groups=[]
    for axis in(0,1):
        values=sorted({F(p[axis])for p in base['pts']})
        groups.extend({'axis':axis,'coordinate':str(v),'members':[j for j,p in enumerate(base['pts'])if F(p[axis])==v]}for v in values)
    assert groups==m['groups']and len(groups)==33 and m['eps']==['1/20000','1/100000','1/500000']
    expected=[(norm(old),{'kind':'baseline'})]
    for quadrant,g,e,sign in product(range(4),groups,map(F,m['eps']),(-1,1)):
        data=norm(old);axis=g['axis'];shift=sign*e
        for j in g['members']:
            k=quadrant*n+j;data['pts'][k][axis]=str((F(data['pts'][k][axis])+shift*P[axis])%(2*P[axis]))
        expected.append((data,{'kind':'quadrant_group_shift','quadrant':quadrant,'group':g,'shift':str(shift)}))
    assert len(expected)==793 and len({digest(d)for d,alias in expected})==793
    assert len(m['tasks'])==len(rows)==len(done)==su['points']==su['completed']==793 and su['errors']==0 and sorted(done,key=lambda x:x['id'])==rows
    assert m['schema']=='quadrant-group-geometry-v1'and m['base_sha256']==sha(bp)and m['target_sha256']==sha(oldp)and m['bits']==40 and m['initial_radius']==5 and m['workers']==6
    assert m['scan_sha256']==sha(t/'computation/scan_quadrant_groups.py')=='bc6d8a19543882b6b84d49cbcc36b51218146bb0b6cc2c70d575cd765cc6d61c'and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')
    assert m['source_prior_manifest_sha256']==sha(prevraw/'manifest.json')and m['source_prior_results_sha256']==sha(prevraw/'results.json')and digest(m['tasks'])==m['task_list_sha256']==su['manifest_task_list_sha256']and sha(raw/'results.json')==su['results_sha256']
    pm=read(prevraw/'manifest.json');pr=read(prevraw/'results.json');prior={x['input_sha256']:(x,z)for x in pm['tasks']for z in pr if x['id']==z['id']};byid={z['id']:z for z in rows};tasks={z['id']:z for z in m['tasks']};assert len(byid)==len(tasks)==793
    actual=Counter();reused=[];kappa=max(w/F(base['pts'][j%n][2])for j,w in enumerate(weights(old)))
    for ix,(task,(data,alias))in enumerate(zip(m['tasks'],expected)):
        assert task['id']==f'{ix:04d}'and norm(task['input'])==data and task['aliases']==[alias]and digest(task['input'])==task['input_sha256']and weights(data)==weights(old)
        assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in data['pts']);delta=displacement(base,data);row=byid[task['id']];R=F(row['tail_radius']);bound=kappa*sum(weights(base))/(P[0]*P[1])*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
        assert 'error'not in row and row['input_sha256']==task['input_sha256']and delta==F(task['delta'])==F(row['delta'])and R==5 and kappa==F(row['kappa'])and bound==F(row['tail_bound'])<1 and F(row['inner_ratio'])<=F(row['outer_ratio'])
        assert row.get('reused',False)==(task['input_sha256']in prior)
        if row.get('reused',False):
            pt,pz=prior[task['input_sha256']];assert norm(pt['input'])==data and row['reused_from_id']==pt['id']=='0110'and row['reused_results_sha256']==sha(prevraw/'results.json')
            assert {k:v for k,v in row.items()if k not in('id','reused','reused_from_id','reused_results_sha256')}=={k:v for k,v in pz.items()if k!='id'};reused.append(task['id'])
        else:
            actual[task['input_sha256']]+=1;cc=row['box_counts'];assert cc['boxes']==cc['zero_mass']+cc['span_rejected']+cc['root_checked']and row['inner_rectangles_unique']<=cc['root_checked']-cc['inner_empty']and row['outer_rectangles_unique']<=cc['root_checked']-cc['outer_empty']
    assert reused==['0000']and len(actual)==sum(actual.values())==su['new_calls']==792 and su['reused']==1
    win=max(rows,key=lambda z:F(z['inner_ratio']));assert win['id']==su['top'][0]['id']=='0360'and norm(tasks['0360']['input'])==norm(final)and sha(raw/'candidate_0360.json')==sha(ip)
    # The extra serial baseline evaluation supplies diagnostics; it is not a new search row.
    baseline=diag['baseline'];assert baseline['id']=='0000'and baseline['input_sha256']==digest(old)
    for key in('delta','kappa','tail_radius','tail_bound','inner_ratio','outer_ratio','window_x_coordinates','window_y_coordinates'):assert baseline[key]==byid['0000'][key]
    differences=[]
    for z in rows:
        if 'box_counts'not in z:continue
        differences.append({'id':z['id'],'box_count_changes':{k:z['box_counts'][k]-baseline['box_counts'][k]for k in z['box_counts']},'inner_rectangles_change':z['inner_rectangles_unique']-baseline['inner_rectangles_unique'],'outer_rectangles_change':z['outer_rectangles_unique']-baseline['outer_rectangles_unique'],'window_x_change':z['window_x_coordinates']-baseline['window_x_coordinates'],'window_y_change':z['window_y_coordinates']-baseline['window_y_coordinates']})
    assert differences==diag['changes']and len(differences)==792
    # Original complete-period kernel certificates, not a geometry rerun in this review.
    p,ind=map(read,(r/'primary_old_R10.json',r/'independent_old_R10.json'));assert sha(raw/'primary_old_R10.json')==sha(r/'primary_old_R10.json')and sha(raw/'independent_old_R10.json')==sha(r/'independent_old_R10.json')and p['input_sha256']==ind['input_sha256']==sha(ip)
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
    lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));assert lo==F(ind['weak_ratio_lower']['fraction'])==F(win['inner_ratio'])and up==F(ind['weak_ratio_upper']['fraction'])==F(win['outer_ratio'])and lo<=up
    assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(ind['area_upper']['fraction'])==up*mass
    assert p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==10
    assert [F(p[k])for k in('period_x','period_y','level','mass')]==[px,py,lam,mass]and list(map(F,ind['periods']))==[px,py]and F(ind['lambda'])==lam and F(ind['base_mass'])==mass
    tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])<1 and q==3
    assert p['counts']['boxes']==ind['counts']['all_span_admissible_boxes']==44361208 and p['window_atoms_distinct']==ind['lifted_atoms']==3953 and p['inner_rectangles_unique']==ind['counts']['unique_inner_rectangles']==24673 and p['outer_rectangles_unique']==ind['counts']['unique_outer_rectangles']==24673
    prevup=F(read(t/'computation/results/round025/primary_final.json')['outer_ratio']);gain=lo-prevup;assert gain>0
    # R24 is an independently audited upper bound for one restricted family only.
    familyfile=t/'computation/results/round024/branch_tree.json';assert sha(familyfile)=='90b6ad00bff70c4346a91063f7e89588bde8e6e6e1c95e13a9b6699a69781108';family=read(familyfile);familyupper=F(family['upper']);familygain=lo-familyupper;assert familyupper==F('5787126969456225479101126541636557987858604811893699578374689/3398607813753623035386245175756010289388060672000000000000000')and familygain>0
    assert familygain==F('570672432886948629541926867364550889447667999850861/16993039068768115176931225878780051446940303360000000000000')
    for n in('manifest','results'):assert gzip.decompress((r/(n+'.json.gz')).read_bytes())==(raw/(n+'.json')).read_bytes()
    assert read(r/'summary.json')==su and read(r/'geometry_diagnostics.json')==diag and sha(r/'completed.jsonl')==sha(raw/'completed.jsonl')
    wall=Decimal(str(su['wall_seconds']));cpu=Decimal(str(su['children_user_cpu_seconds']))+Decimal(str(su['children_system_cpu_seconds']));assert f'{wall:.3f}'=='680.948'and f'{cpu:.3f}'=='3385.655'
    proof=t/'proofs/round026-lower.md';report=t/'rounds/026.md';prooftext=proof.read_text();assert str(lo.numerator)in prooftext and str(lo.denominator)in prooftext and str(up.numerator)in prooftext and str(up.denominator)in prooftext
    assert str(tail.numerator)in prooftext and str(tail.denominator)in prooftext and str(familygain.numerator)in prooftext and str(familygain.denominator)in prooftext
    reporttext=report.read_text();assert all(sha(z)in reporttext for z in(raw/'manifest.json',raw/'results.json',r/'manifest.json.gz',r/'results.json.gz'))
    familyreview=t/'reviews/affine-mass-envelope-review.md';familyaudit=t/'computation/results/round024/final-checkpoint-check.json'
    assert sha(familyreview)=='e895cc868858e9e25a3f3a0118182c465ff711c2ca848da5345fa471b5928c49'and sha(familyaudit)=='f1ac339db17db3e2c8aadc4a6c020a1bed67f3141ec44dd4e703eb8d8d6c13b5'
    assert sha(proof)=='3a59abe184c7cc70108ea73d934fca71d1295d011046d1912ed0a729977574e1'and sha(report)=='f0ee9ebfd3215bc9c57504a635be69fe60d933a0cd040dbd4eaf31ba39e52751'
    assert sha(familyreview)in prooftext and sha(familyaudit)in prooftext
    files=[ip,proof,report,familyreview,familyaudit,t/'computation/build_round026_ledger.py',familyfile,t/'computation/scan_quadrant_groups.py',*sorted(r.glob('*'))]
    out={'status':'passed','input_sha256':sha(ip),'canonical_input_sha256':digest(final),'mass':str(mass),'positive_masses':204,'changes_zero_based':changes,'coordinate_group_count':len(groups),'winner_quadrant':1,'winner_members':[28,29],'winner_axis':'x','winner_shift_fraction_of_R8_period':'1/500000','tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R25_full_upper':str(gain),'R24_restricted_family_upper':str(familyupper),'R24_certificate_sha256':sha(familyfile),'lower_minus_R24_restricted_family_upper':str(familygain),'R24_scope':'Fixed R14 positions, D=[-1,1] x [0,10] affine mass family only; does not upper-bound general C2','ledger_counts':{'rows':793,'unique_inputs':793,'new_search_evaluations':792,'reuse':1,'extra_serial_baseline_diagnostic_evaluations':1},'diagnostic_recorded_seconds':baseline['elapsed_seconds'],'diagnostic_delta_rows':len(differences),'all_inputs_independently_reconstructed':True,'compressed_archives_match_raw':True,'geometry_or_search_rerun':False,'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'wall_seconds':str(wall),'child_cpu_seconds':str(cpu),'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'audit_helper_sha256':sha(Path(__file__).with_name('review_round017_ledger.py')),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    # Bind producer integration/ledger after all independent reconstruction.
    assert c['input_sha256']==sha(ip)and F(c['lower'])==lo and F(c['upper'])==up and F(c['tail_bound'])==tail and c['q']==q
    assert c['primary_sha256']==sha(r/'primary_old_R10.json')and c['independent_sha256']==sha(r/'independent_old_R10.json')and c['search_rows']==793 and c['search_new_calls']==792 and c['search_reused']==1 and c['additional_diagnostic_baseline_calls']==1 and c['additional_diagnostic_seconds']==baseline['elapsed_seconds']and c['changed_label_indices_zero_based']==[79,80]and F(c['gain_over_r25_upper'])==gain and F(c['gain_over_r24_family_upper'])==familygain
    assert c['wall_seconds']==su['wall_seconds']and c['child_cpu_seconds']==su['children_user_cpu_seconds']+su['children_system_cpu_seconds']
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
