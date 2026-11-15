"""Independent Fraction reconstruction of the complete fixed-mass R25 grid."""
from review_round017_ledger import *
import gzip

def labelled_positions(base,masses,A,B):
    # Independent formula, no production geometry-builder import.
    P=[F(base['Px']),F(base['Py'])];pts=[]
    for cx,cy in product((0,1),repeat=2):
        for j,row in enumerate(base['pts']):
            shift=F(0)
            if cy==1 and 27<=j<=49 and j%2==1:shift+=A*3*P[1]/400000
            if cx==cy==1 and j==50:shift+=B*P[1]/400000
            pts.append([str((F(row[0])+cx*P[0])%(2*P[0])),str((F(row[1])+cy*P[1]+shift)%(2*P[1])),str(masses[len(pts)])])
    return {'Px':str(2*P[0]),'Py':str(2*P[1]),'lam':str(F(base['lam'])),'pts':pts}

def displacement(base,target):
    P=[F(base['Px']),F(base['Py'])];delta=F(0);n=len(base['pts'])
    for k,row in enumerate(target['pts']):
        cx,cy=divmod(k//n,2);j=k%n
        for axis,copy in ((0,cx),(1,cy)):
            d=(F(row[axis])-F(base['pts'][j][axis])-copy*P[axis])%(2*P[axis]);delta=max(delta,min(d,2*P[axis]-d))
    return delta

def certificate(t,r,ip,final,c,previous_upper):
    p,ind=map(read,(r/'primary_final.json',r/'independent_final.json'));px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));R=F(10);tail=mass/(px*py)*(1+px/(2*R))*(1+py/(2*R));q=math.ceil(R/min(px,py))+1;lo,up=map(F,(p['inner_ratio'],p['outer_ratio']))
    assert sha(ip)==p['input_sha256']==ind['input_sha256']==c['input_sha256']
    assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'
    assert ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
    assert p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==R and q==c['q']==3
    assert [F(p[k])for k in('period_x','period_y','level','mass')]==[px,py,lam,mass] and list(map(F,ind['periods']))==[px,py] and F(ind['lambda'])==lam and F(ind['base_mass'])==mass
    assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])==F(c['tail_bound'])<lam==1
    assert lo==F(ind['weak_ratio_lower']['fraction'])==F(c['lower']) and up==F(ind['weak_ratio_upper']['fraction'])==F(c['upper']) and lo<=up
    assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==mass*lo/lam and F(p['outer_area'])==F(ind['area_upper']['fraction'])==mass*up/lam
    assert lo-previous_upper==F(c['gain_over_r23_upper'])>0
    assert p['input_atoms']==ind['base_atoms']==len(final['pts'])==204 and p['window_atoms_distinct']==ind['lifted_atoms'] and p['counts']['boxes']==ind['counts']['all_span_admissible_boxes'] and p['inner_rectangles_unique']==ind['counts']['unique_inner_rectangles'] and p['outer_rectangles_unique']==ind['counts']['unique_outer_rectangles']
    return p,ind,lo,up,tail,q

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round025'
    ip=t/'computation/inputs/r025_final_0110.json';bp=t/'computation/inputs/r008_final_0611.json';oldp=t/'computation/inputs/r023_final_1051.json';final,base,old=map(read,(ip,bp,oldp));c=read(r/'integration_checks.json');ledger=read(r/'search_ledger.json');m=read(raw/'manifest.json');rows=read(raw/'results.json');su=read(raw/'summary.json');done=[json.loads(s)for s in(raw/'completed.jsonl').read_text().splitlines()]
    assert sha(ip)=='34a9569786741ddf2a61cc52fd85725bf8efdccf46306e34c1d8d340be49cd0c' and sha(oldp)==c['baseline_r23_sha256']=='4fb1461289fd2b1ef52ac7d30c3d5cd0bd5b72c102f7ccdbeb01be6706baf0ee'
    assert weights(final)==weights(old) and len(base['pts'])==51 and len(final['pts'])==204 and min(weights(final))>0
    mass=sum(weights(final));P=[F(base['Px']),F(base['Py'])];px,py=F(final['Px']),F(final['Py']);assert mass==4*sum(weights(base))==F(c['mass']) and (px,py)==(2*P[0],2*P[1]) and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts'])
    assert norm(final)==labelled_positions(base,weights(old),F(9,8),F(27,64)) and list(map(F,c['parameters']))==[F(9,8),F(27,64)]
    changes={name:[j for j,(u,v)in enumerate(zip(old['pts'],final['pts']))if F(u[k])!=F(v[k])]for k,name in enumerate(('x','y','weight'))};assert changes==c['changes_zero_based']
    oldup=F(read(t/'computation/results/round023/primary_final.json')['outer_ratio']);p,ind,lo,up,tail,q=certificate(t,r,ip,final,c,oldup)
    assert sha(raw/'primary_old_R10.json')==sha(r/'primary_final.json') and sha(raw/'candidate_0110.json')==sha(ip)
    assert m['schema']=='fixed-mass-two-position-modes-v1'and m['base_sha256']==sha(bp)and m['mass_input_sha256']==sha(oldp)and m['bits']==40 and m['initial_tail_radius']==5 and m['workers']==6
    assert m['scan_sha256']==ledger['scan_sha256']==sha(t/'computation/scan_fixedmass_positions.py')and m['geometry_builder_sha256']==sha(t/'computation/scan_supercell_joint.py')and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')
    assert sha(raw/'manifest.json')==ledger['manifest_sha256']and sha(raw/'results.json')==ledger['results_sha256']==su['results_sha256']and digest(m['tasks'])==m['task_list_sha256']==ledger['task_list_sha256']==su['manifest_task_list_sha256']
    assert len(m['tasks'])==len(rows)==len(done)==su['points']==su['completed']==ledger['points']==c['result_rows']==289 and su['errors']==0 and sorted(done,key=lambda x:x['id'])==rows
    assert read(r/'summary.json')==su
    for n in('manifest','results'):assert gzip.decompress((r/(n+'.json.gz')).read_bytes())==(raw/(n+'.json')).read_bytes()
    tasks={x['id']:x for x in m['tasks']};byid={x['id']:x for x in rows};assert len(tasks)==len(byid)==289
    combos=set();hashes=set()
    for tid,task in tasks.items():
        data=task['input'];assert digest(data)==task['input_sha256'] and task['input_sha256']not in hashes;hashes.add(task['input_sha256']);assert weights(data)==weights(old)
        assert len(task['aliases'])==1
        for alias in task['aliases']:
            i,j=alias['i'],alias['j'];assert (i,j)not in combos;combos.add((i,j));A,B=F(19,16)+F(i,32),F(27,64)+F(j,32);assert [F(alias[k])for k in('a','b')]==[A,B]and norm(data)==labelled_positions(base,weights(old),A,B)
        delta=displacement(base,data);row=byid[tid];kappa=max(w/F(base['pts'][j%51][2])for j,w in enumerate(weights(data)));R=F(row['tail_radius']);bound=kappa*sum(weights(base))/(P[0]*P[1])*(1+(P[0]+2*delta)/(2*R))*(1+(P[1]+2*delta)/(2*R))
        assert 'error'not in row and not row.get('reused',False)and row['input_sha256']==task['input_sha256']and delta==F(task['delta'])==F(row['delta'])and R==5 and kappa==F(row['kappa']) and bound==F(row['tail_bound'])<1 and F(row['inner_ratio'])<=F(row['outer_ratio'])
    assert combos==set(product(range(-8,9),repeat=2))and ledger['new_calls']==c['worker_evaluations']==289 and ledger['reused']==c['reused_records']==0
    win=max(rows,key=lambda x:F(x['inner_ratio']));assert win['id']==ledger['winner']['id']==su['top'][0]['id']=='0110'and norm(tasks[win['id']]['input'])==norm(final)and F(win['inner_ratio'])==lo and F(win['outer_ratio'])==up
    for k,v in win.items():assert ledger['winner'][k]==su['top'][0][k]==v
    assert ledger['winner']['aliases']==tasks[win['id']]['aliases']
    wall=Decimal(str(su['wall_seconds']));cpu=Decimal(str(su['children_user_cpu_seconds']))+Decimal(str(su['children_system_cpu_seconds']));assert f'{wall:.3f}'=='171.429'and f'{cpu:.3f}'=='999.249'
    proof=t/'proofs/round025-lower.md';report=t/'rounds/025.md';text=proof.read_text();fr=[F(int(u),int(v))for u,v in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',text)];assert fr==[lo,up]and str(tail)in text and str(lo-oldup)in text and sha(ip)in text and digest(final)in text
    files=[ip,proof,report,t/'computation/scan_fixedmass_positions.py',t/'computation/scan_supercell_joint.py',*sorted(r.glob('*'))]
    out={'status':'passed','input_sha256':sha(ip),'canonical_input_sha256':digest(final),'mass':str(mass),'positive_masses':204,'parameters':['9/8','27/64'],'changes_zero_based':changes,'tail_bound':str(tail),'q':q,'lower':str(lo),'upper':str(up),'lower_minus_R23_full_upper':str(lo-oldup),'ledger_counts':{'rows':289,'unique_inputs':len(hashes),'actual':289,'reuse':0},'all_inputs_independently_reconstructed':True,'compressed_archives_match_raw':True,'geometry_or_search_rerun':False,'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'wall_seconds':str(wall),'child_cpu_seconds':str(cpu),'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'audit_helper_sha256':sha(Path(__file__).with_name('review_round017_ledger.py')),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
