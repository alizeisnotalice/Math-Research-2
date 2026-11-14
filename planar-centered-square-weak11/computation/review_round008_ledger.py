"""Read-only R8 exact proof/ledger check from the fixed R6 base.
The stale coordinate value in a direction alias is intentionally not used.
No geometry calculator or search is run.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from datetime import datetime,timezone
import argparse,hashlib,json,math,platform,re,sys,time

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def norm(d):return {'Px':str(F(d['Px'])),'Py':str(F(d['Py'])),'lam':str(F(d['lam'])),
                    'pts':[[str(F(v)) for v in p] for p in d['pts']]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True)
    ap.add_argument('--raw-search',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.perf_counter();t=args.topic
    ip=t/'computation/inputs/r008_final_0611.json';bp=t/'computation/inputs/r006_final_0015.json'
    d,b=read(ip),read(bp);base=norm(b);res=t/'computation/results/round008';r6=t/'computation/results/round006'
    p,i,c=map(read,(res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json'))
    old=read(r6/'primary_final.json');ledger=read(res/'search_ledger.json')
    assert len(ledger)==3 and len(d['pts'])==len(b['pts'])==51
    assert sha(ip)==p['input_sha256']==i['input_sha256']==c['input_sha256']
    assert sha(bp)==c['baseline_input_sha256']
    px,py,lam=map(F,(d['Px'],d['Py'],d['lam']));mass=sum(F(q[2]) for q in d['pts'])
    assert all(F(d[k])==F(b[k]) for k in ('Px','Py','lam')) and lam==1
    assert mass==sum(F(q[2]) for q in b['pts']) and all(F(q[2])>0 for q in d['pts'])
    assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
    tail=mass/(px*py)*(1+px/10)*(1+py/10);q=math.ceil(F(5)/min(px,py))+1
    assert tail==F(c['tail_bound'])==F(p['tail_bound'])==F(i['tail_bound_at_cap']['fraction'])<1 and q==c['q']==3
    lo,up=F(p['inner_ratio']),F(p['outer_ratio']);margin=lo-F(old['outer_ratio'])
    assert lo==F(i['weak_ratio_lower']['fraction'])==F(c['lower'])
    assert up==F(i['weak_ratio_upper']['fraction'])==F(c['upper'])
    assert F(p['inner_area'])==F(i['area_lower']['fraction'])==lo*mass
    assert F(p['outer_area'])==F(i['area_upper']['fraction'])==up*mass
    assert margin==F(c['difference_over_baseline_outer'])>0
    assert p['bits']==i['bits']==40 and F(c['epsilon'])==F(1,2**41)
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')
    assert i['script_sha256']==sha(t/'computation/independent_periodic_halo.py')
    dirs=ledger[0]['directions']
    assert dirs[0]['kind']=='mass' and dirs[0]['group']['members']==[0] and dirs[0]['anchor']==1
    assert dirs[0]['sign']==-1 and F(dirs[0]['eps'])==F(1,100000)
    assert dirs[1]['kind']=='coordinate' and dirs[1]['group']['axis']==0 and dirs[1]['group']['members']==[2,3]
    assert dirs[1]['sign']==1 and F(dirs[1]['eps'])==F(1,100000)
    source_summaries=[r6/'epoch7_summary.json',r6/'epoch5_summary.json']
    assert [sha(z) for z in source_summaries]==ledger[0]['direction_summaries_sha256']
    m0,m1=F(base['pts'][0][2]),F(base['pts'][1][2]);assert min(m0,m1)==m0
    assert F(dirs[1]['group']['coordinate'])!=F(base['pts'][2][0])
    def rebuild(coeffs):
        a,z=map(F,coeffs);out={k:v for k,v in base.items() if k!='pts'}
        pts=[p.copy() for p in base['pts']]
        amount=a*m0/100000
        pts[0][2]=str(m0-amount);pts[1][2]=str(m1+amount)
        for j in (2,3):pts[j][0]=str((F(base['pts'][j][0])+z*px/100000)%px)
        out['pts']=pts
        assert all(F(p[2])>0 for p in pts)
        assert sum(F(p[2]) for p in pts)==mass
        return out
    expected_coeffs=[F(-1,32),F(69,64)]
    assert rebuild(expected_coeffs)==norm(d) and list(map(F,c['winner_coefficients']))==expected_coeffs
    counts={'x':[],'y':[],'weight':[]};changes=[]
    for index,(a,z) in enumerate(zip(b['pts'],d['pts'])):
        for axis,name in enumerate(('x','y','weight')):
            if F(a[axis])!=F(z[axis]):
                counts[name].append(index);changes.append({'zero_based_atom':index,'field':name,
                  'before':str(F(a[axis])),'after':str(F(z[axis])),'delta':str(F(z[axis])-F(a[axis]))})
    assert counts==c['coordinate_changes_indices_zero_based']=={'x':[2,3],'y':[],'weight':[0,1]}
    wall=Decimal(0);cpu=Decimal(0);total=0;reconstructed=0;rows=[]
    expected_centers=[[F(0),F(0)],[F(0),F(1)],[F(-1,16),F(17,16)]]
    expected_steps=[F(1,4),F(1,16),F(1,64)]
    for batch,row in enumerate(ledger):
        raw=args.raw_search/f'batch{batch}';manifest=read(raw/'manifest.json')
        results=read(raw/'results.json');summary=read(res/f'batch{batch}_summary.json')
        completed=[json.loads(line) for line in (raw/'completed.jsonl').read_text().splitlines()]
        assert row['input_sha256']==manifest['input_sha256']==sha(bp)
        assert row['directions']==manifest['directions']==dirs
        assert row['points']==row['completed']==summary['points']==summary['completed']==1089
        assert row['errors']==summary['errors']==row['rejected']==summary['rejected']==0
        assert row['workers']==12 and row['bits']==40
        assert sha(raw/'manifest.json')==row['manifest_sha256']
        assert sha(raw/'results.json')==row['results_sha256']==summary['results_sha256']
        assert digest(manifest['tasks'])==row['task_list_sha256']==summary['manifest_task_list_sha256']
        assert len(results)==len(completed)==len(manifest['tasks'])==1089
        assert len({z['id'] for z in results})==1089
        assert sorted(completed,key=lambda z:z['id'])==sorted(results,key=lambda z:z['id'])
        center=list(map(F,manifest['center']));step=F(manifest['step'])
        assert center==list(map(F,row['center']))==expected_centers[batch] and step==F(row['step'])==expected_steps[batch]
        ij=set();tasks={z['id']:z for z in manifest['tasks']};center_id=None
        for task in tasks.values():
            assert digest(task['input'])==task['input_sha256']
            for alias in task['aliases']:
                aa,bb=alias['i'],alias['j'];assert -16<=aa<=16 and -16<=bb<=16;ij.add((aa,bb))
                coeff=[center[0]+aa*step,center[1]+bb*step]
                assert list(map(F,alias['coefficients']))==coeff
                assert rebuild(coeff)==norm(task['input'])
                reconstructed+=1
                if aa==bb==0:center_id=task['id']
        assert len(ij)==1089
        assert all('error' not in z and z['input_sha256']==tasks[z['id']]['input_sha256']
                   and F(z['inner_ratio'])<=F(z['outer_ratio']) and z['tail_radius']==5 for z in results)
        winner=max(results,key=lambda z:F(z['inner_ratio']));center_result=next(z for z in results if z['id']==center_id)
        assert winner['id']==row['winner_id']==summary['top'][0]['id']
        assert F(winner['inner_ratio'])==F(row['winner_inner_ratio'])>F(center_result['outer_ratio'])==F(row['center_outer_ratio'])
        assert tasks[winner['id']]['input_sha256']==row['winner_input_canonical_sha256']
        assert tasks[winner['id']]['aliases']==row['winner_aliases']
        winner_coeff=list(map(F,row['winner_aliases'][0]['coefficients']))
        if batch<2:assert winner_coeff==expected_centers[batch+1]
        else:assert winner_coeff==expected_coeffs and norm(tasks[winner['id']]['input'])==norm(d)
        wall+=Decimal(str(row['wall_seconds']));cpu+=Decimal(str(row['children_user_cpu_seconds']))+Decimal(str(row['children_system_cpu_seconds']))
        total+=len(results);rows.append({'batch':batch,'tasks_reconstructed':len(tasks),'winner_coefficients':list(map(str,winner_coeff)),
                                       'strict_gain_over_center':str(F(winner['inner_ratio'])-F(center_result['outer_ratio']))})
    assert total==reconstructed==3267 and f'{wall:.6f}'=='511.075172' and f'{cpu:.6f}'=='5630.669179'
    frac=lambda s:[F(int(a.strip()),int(b.strip())) for a,b in re.findall(r'\\frac\{([^{}]+)\}\s*\{([^{}]+)\}',s)]
    proof=t/'proofs/round008-lower.md';report=t/'rounds/008.md'
    assert frac(proof.read_text())==[mass,tail,lo,up,margin] and frac(report.read_text())==[lo]
    r1=Path('/Volumes/data/CodexProjects/T2/planar-centered-square-weak11/proofs/periodic-complete-halo.md')
    assert sha(r1)=='7571ad4fa597d9f152fcfd1249bf07c69d998326abf30c44a9abaca1fe18987e'
    files=[ip,proof,report,res/'primary_final.json',res/'independent_final.json',res/'integration_checks.json',res/'search_ledger.json']
    out={'status':'passed','scope':'fixed-base transformations and saved arithmetic/ledger; no geometry or search rerun',
         'input_sha256':sha(ip),'fixed_base_R6_sha256':sha(bp),'winner_coefficients':list(map(str,expected_coeffs)),
         'stale_coordinate_label':dirs[1]['group']['coordinate'],'actual_R6_x2':base['pts'][2][0],
         'coordinate_label_ignored_in_transform':True,'mass_transfer_from_1_to_0':str(m0/F(3200000)),
         'x2_x3_displacement':str(F(69,6400000)*px),'positive_masses':51,'mass':str(mass),
         'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R6_full_upper':str(margin),
         'exact_input_changes':changes,'all_manifest_inputs_reconstructed':reconstructed,'completion_records_checked':total,
         'batch_checks':rows,'wall_seconds_sum':str(wall),'child_cpu_seconds_sum':str(cpu),
         'files_sha256':{str(p):sha(p) for p in files},'audit_script_sha256':sha(Path(__file__)),
         'verified_utc':datetime.now(timezone.utc).isoformat(),'python':sys.version,'platform':platform.platform(),
         'elapsed_seconds':time.perf_counter()-start,'command':[sys.executable,*sys.argv]}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('exact_input_changes','files_sha256')},indent=2))

if __name__=='__main__':main()
