"""Independent Fraction reconstruction of R41 selection, all inputs and certificate ledger.
No periodic geometry/kernel evaluation is performed by this review.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from itertools import combinations
from datetime import datetime,timezone
import argparse,gzip,hashlib,json,math,os,re,sys,time

def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def norm(x): return {**{k:str(F(x[k]))for k in ('Px','Py','lam')},'pts':[[str(F(v))for v in q]for q in x['pts']]}
def weights(x): return [F(q[2])for q in x['pts']]
def delta(base,data):
 P=[F(base[k])for k in ('Px','Py')]; n=len(base['pts']); ans=F(0)
 for i,row in enumerate(data['pts']):
  cx,cy=divmod(i//n,2)
  for axis,copy in ((0,cx),(1,cy)):
   d=(F(row[axis])-F(base['pts'][i%n][axis])-copy*P[axis])%(2*P[axis]);ans=max(ans,min(d,2*P[axis]-d))
 return ans

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round041';prior=raw.parent.parent/'round026/quadrant-groups'
 assert os.path.ismount('/Volumes/data')
 bp=t/'computation/inputs/r008_final_0611.json';op=t/'computation/inputs/r035_final_0134.json';ip=t/'computation/inputs/r041_final_0036.json';base,old,final=map(read,(bp,op,ip));m=read(raw/'manifest.json');rows=read(raw/'results.json');su=read(raw/'summary.json');done=[json.loads(s)for s in(raw/'completed.jsonl').read_text().splitlines()];c=read(r/'integration_checks.json')
 assert sha(ip)=='ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3' and sha(op)=='5021c54704af45b3ef8289c67a245926473434640372a934b2f15c4485f0c026'
 n=len(base['pts']);P=[F(base[k])for k in ('Px','Py')];px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));assert n==51 and len(old['pts'])==len(final['pts'])==204 and min(weights(final))>0 and weights(final)==weights(old)and mass==4*sum(weights(base))and (px,py)==tuple(2*v for v in P)and lam==1
 assert all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts'])
 assert all(F(final[k])==F(old[k])for k in ('Px','Py','lam'))
 changes={name:[i for i,(u,v)in enumerate(zip(old['pts'],final['pts']))if F(u[axis])!=F(v[axis])]for axis,name in enumerate(('x','y','mass'))};assert changes=={'x':[77,78],'y':[],'mass':[]}
 assert all(F(final['pts'][i][0])-F(old['pts'][i][0])==F(-328734960629,120000000000000000) for i in changes['x'])
 assert m['schema']=='r41-cross-group-probe-v1' and m['base_sha256']==sha(bp)and m['target_sha256']==sha(op)and m['bits']==40 and m['radius']==8 and m['workers']==12
 assert m['script_sha256']==sha(t/'computation/scan_cross_group_r41.py')=='2beda1ad4769ab459ebc319ded387c00ecc46e7b78b27dea3f6fedc840374a54'
 assert m['kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')=='a697093e30172da71528dc3668add00ea4083811cc2b69ef3a0fc9942eded5e3'
 assert m['source_manifest_sha256']==sha(prior/'manifest.json')=='4445cfc88d7ba633bd8b3d495a1aa90a119d542a751b23fe917931e85953f691' and m['source_results_sha256']==sha(prior/'results.json')=='0e7633d75509ef95a353a05901d45b8a3605c4199e725d28e93d6974fd7759fb'
 selectionold=read(t/'computation/inputs/r025_final_0110.json');pm=read(prior/'manifest.json');pr=read(prior/'results.json');pt={z['id']:z for z in pm['tasks']};pz={z['id']:z for z in pr};assert pz['0000']['input_sha256']==digest(selectionold)
 # Exact ranking and first occurrence tie rule, followed by quadrant/axis/member/sign deduplication.
 ranked=sorted((z for z in pr if F(z['inner_ratio'])>F(pz['0000']['outer_ratio'])),key=lambda z:F(z['inner_ratio']),reverse=True);directions=[];seen=set()
 for z in ranked:
  src=pt[z['id']];alias=src['aliases'][0];g=alias['group'];axis=g['axis'];members=[j for j,p in enumerate(base['pts'])if F(p[axis])==F(g['coordinate'])];assert members==g['members'] and src['input_sha256']==digest(src['input'])==z['input_sha256']
  expected=norm(selectionold)
  for j in members:expected['pts'][alias['quadrant']*n+j][axis]=str((F(expected['pts'][alias['quadrant']*n+j][axis])+F(alias['shift'])*P[axis])%(2*P[axis]))
  assert expected==norm(src['input'])
  key=(alias['quadrant'],axis,tuple(members),F(alias['shift'])>0)
  if key not in seen:seen.add(key);directions.append({'source_id':z['id'],'alias':alias})
 assert len(ranked)==m['strict_improving_count']==61 and len(directions)==m['distinct_direction_count']==45
 dirs=[{'quadrants':[q],'members':g}for g in ([26,27],[32,33])for q in range(4)];assert m['directions']==dirs
 expected={}
 for direction in dirs:
  for j in range(-12,13):
   d=norm(old)
   for q in direction['quadrants']:
    for label in direction['members']:
     k=q*n+label;d['pts'][k][0]=str((F(d['pts'][k][0])+F(j,1200000)*P[0])%px)
   h=digest(d)
   if h not in expected:expected[h]={'id':'','input_sha256':h,'input':d,'delta':str(delta(base,d)),'combinations':[]}
   expected[h]['combinations'].append({'quadrants':direction['quadrants'],'members':direction['members'],'step':str(F(j,1200000))})
 h=digest(old);baseline=expected.pop(h);expected={h:baseline,**expected}
 for j,item in enumerate(expected.values()):item['id']=f'{j:04d}'
 assert len(expected)==193 and sum(len(x['combinations'])for x in expected.values())==200
 assert len(rows)==len(m['tasks'])==len(done)==su['points']==su['completed']==193 and su['new_calls']==192 and su['reused']==1 and su['errors']==0
 assert sorted(done,key=lambda x:x['id'])==rows and digest(m['tasks'])==m['task_list_sha256']and sha(raw/'results.json')==su['results_sha256']
 assert len({z['id']for z in rows})==len({z['input_sha256']for z in rows})==193
 kappa=max(w/F(base['pts'][i%n][2])for i,w in enumerate(weights(old)));assert kappa==F(25021,25000)
 byid={z['id']:z for z in rows};oldmatch=[]
 for task,ex in zip(m['tasks'],expected.values()):
  assert task==ex
  d=task['input'];row=byid[task['id']];dd=delta(base,d);R=F(row['tail_radius']);bound=kappa*sum(weights(base))/(P[0]*P[1])*(1+(P[0]+2*dd)/(2*R))*(1+(P[1]+2*dd)/(2*R))
  assert weights(d)==weights(old)and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in d['pts'])
  assert 'error'not in row and row.get('reused',False)==(task['id']=='0000')and row['input_sha256']==task['input_sha256']and F(row['delta'])==F(task['delta'])==dd and R==8 and F(row['kappa'])==kappa and F(row['tail_bound'])==bound<lam and F(row['inner_ratio'])<=F(row['outer_ratio'])and row['elapsed_seconds']>0
  cc=row['box_counts'];assert cc['boxes']==cc['zero_mass']+cc['span_rejected']+cc['root_checked']and 0<=cc['inner_empty']<=cc['root_checked']and 0<=cc['outer_empty']<=cc['root_checked']
  if task['input_sha256'] in {z['input_sha256']for z in pr}:oldmatch.append(task['id'])
 # The single reused row is the entire old summary winner, including its old combination label.
 prevraw=raw.parent.parent/'round035/axis-wall-refine';prevsummary=read(prevraw/'summary.json');prevrows=read(prevraw/'results.json');previouswinner=prevsummary['top'][0]
 assert m['source_r35_summary_sha256']==sha(prevraw/'summary.json')and previouswinner['id']=='0134'
 assert byid['0000']=={**previouswinner,'id':'0000','reused':True,'source_results_sha256':sha(prevraw/'results.json')}
 prevwinnerrow=next(z for z in prevrows if z['id']==previouswinner['id'])
 assert {k:v for k,v in previouswinner.items()if k!='combinations'}==prevwinnerrow and previouswinner['input_sha256']==digest(old)
 assert m['tasks'][0]['combinations']==[{'quadrants':d['quadrants'],'members':d['members'],'step':'0'}for d in dirs] and byid['0000']['combinations']==[{'quadrants':[2],'step':'9/2500000'}]
 assert sum(not z.get('reused',False)for z in rows)==192
 ordered=sorted(rows,key=lambda z:F(z['inner_ratio']),reverse=True)
 for row,top in zip(ordered[:8],su['top']):assert top=={**row,'combinations':m['tasks'][int(row['id'])]['combinations']}
 win=ordered[0];assert win['id']=='0036' and m['tasks'][36]['combinations']==[{'quadrants':[1],'members':[26,27],'step':'-1/1200000'}]and norm(m['tasks'][36]['input'])==norm(final)and sha(raw/'candidate_0036.json')==sha(ip)
 # Independently form all exact slope intervals, including baseline aliases.
 analysis=read(raw/'secant_analysis.json');groups={}
 for task in m['tasks']:
  for alias in task['combinations']:groups.setdefault((tuple(alias['members']),tuple(alias['quadrants'])),[]).append((F(alias['step']),byid[task['id']]))
 expected_analysis={};ndisjoint=0
 for quads,samples in groups.items():
  samples.sort(key=lambda z:z[0]);segments=[]
  for (aa,u),(bb,v)in zip(samples,samples[1:]):
   step=bb-aa;low=(F(v['inner_ratio'])-F(u['outer_ratio']))/step;high=(F(v['outer_ratio'])-F(u['inner_ratio']))/step;assert step>0 and low<=high
   segments.append({'left':str(aa),'right':str(bb),'slope_lower':str(low),'slope_upper':str(high),'slope_lower_float':float(low),'slope_upper_float':float(high)})
  disjoint=[]
  for left,right in zip(segments,segments[1:]):
   if F(left['slope_upper'])<F(right['slope_lower'])or F(right['slope_upper'])<F(left['slope_lower']):disjoint.append({'left_interval':[left['left'],left['right']],'right_interval':[right['left'],right['right']],'left_slope':[left['slope_lower'],left['slope_upper']],'right_slope':[right['slope_lower'],right['slope_upper']]})
  expected_analysis[str(quads)]={'sample_points':len(samples),'segments':segments,'disjoint_adjacent_secants':disjoint};ndisjoint+=len(disjoint)
 assert analysis['directions']==expected_analysis and analysis['manifest_sha256']==sha(raw/'manifest.json')and analysis['results_sha256']==sha(raw/'results.json')and analysis['script_sha256']==sha(t/'computation/analyze_cross_group_r41.py')
 assert sha(r/'secant_analysis.json')==sha(raw/'secant_analysis.json')
 # Match the two already executed complete-period kernels, with their frozen algorithm identities.
 p,ind=map(read,(r/'primary_old_R10.json',r/'independent_old_R10.json'))
 for name in ('primary_old_R10.json','independent_old_R10.json','summary.json','completed.jsonl'):assert sha(raw/name)==sha(r/name)
 assert p['input_sha256']==ind['input_sha256']==sha(ip)
 assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33' and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
 lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));assert lo==F(ind['weak_ratio_lower']['fraction'])==F(win['inner_ratio'])and up==F(ind['weak_ratio_upper']['fraction'])==F(win['outer_ratio'])and lo<=up
 assert F(p['inner_area'])==F(ind['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(ind['area_upper']['fraction'])==up*mass
 assert p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==10
 assert [F(p[k])for k in ('period_x','period_y','level','mass')]==[px,py,lam,mass]and list(map(F,ind['periods']))==[px,py]and F(ind['lambda'])==lam and F(ind['base_mass'])==mass
 tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])<lam and q==3
 assert p['counts']['boxes']==ind['counts']['all_span_admissible_boxes'] and p['window_atoms_distinct']==ind['lifted_atoms'] and p['inner_rectangles_unique']==ind['counts']['unique_inner_rectangles'] and p['outer_rectangles_unique']==ind['counts']['unique_outer_rectangles']
 prevup=F(read(t/'computation/results/round035/primary_old_R10.json')['outer_ratio']);gain=lo-prevup;assert gain>0
 for name in ('manifest','results'):assert gzip.decompress((r/(name+'.json.gz')).read_bytes())==(raw/(name+'.json')).read_bytes()
 proof=t/'proofs/round041-lower.md';report=t/'rounds/041.md';transfer=t/'proofs/periodic-complete-halo.md'
 fr=[F(int(x),int(y))for x,y in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]and sha(proof)=='e84c132a63410709f79f4f6f2b3aef8368e659f0d0a5b80e7e33e2422a4c9f9a'and sha(report)=='6e1bbe0304bea83e8f96093b8fe6b2d33468154f74d19d00ee0c8df4c8efb002'
 assert sha(ip)in proof.read_text()and all(sha(x)in report.read_text()for x in [raw/'manifest.json',raw/'results.json',r/'manifest.json.gz',r/'results.json.gz'])
 assert c['input_sha256']==sha(ip)and c['primary_sha256']==sha(r/'primary_old_R10.json')and c['independent_sha256']==sha(r/'independent_old_R10.json')and c['search_rows']==193 and c['search_new_calls']==192 and c['reused']==1 and F(c['lower'])==lo and F(c['upper'])==up and F(c['gain_over_r35_upper'])==gain and F(c['tail_bound'])==tail and c['q']==3 and c['changed_label_indices_zero_based']==changes['x']
 wall=Decimal(str(su['wall_seconds']));cpu=Decimal(str(su['child_user_cpu_seconds']))+Decimal(str(su['child_system_cpu_seconds']));assert f'{wall:.6f}'=='275.633482'and f'{cpu:.6f}'=='2860.600122'and c['wall_seconds']==su['wall_seconds']and c['child_cpu_seconds']==su['child_user_cpu_seconds']+su['child_system_cpu_seconds']
 files=[bp,op,ip,proof,report,transfer,t/'computation/analyze_cross_group_r41.py',raw/'secant_analysis.json',t/'computation/scan_cross_group_r41.py',t/'computation/build_round041_ledger.py',t/'computation/periodic_halo.py',t/'computation/independent_periodic_halo.py',t/'computation/periodic_halo_dominated.py',t/'computation/scan_quadrant_groups.py',t/'computation/scan_supercell_mass.py',t/'computation/scan_supercell_2x2.py',t/'computation/results/round035/primary_old_R10.json',prior/'manifest.json',prior/'results.json',prevraw/'summary.json',prevraw/'results.json',*sorted(r.glob('*'))]
 out={'status':'passed','scope':'Exact ledger reconstruction and two existing frozen-kernel certificates; no geometry rerun','input_sha256':sha(ip),'canonical_input_sha256':digest(final),'lower':str(lo),'upper':str(up),'lower_minus_R35_full_upper':str(gain),'tail_bound_R10':str(tail),'q':q,'D':str(2**40),'epsilon':str(F(1,2**41)),'finite_L1_factor':'(1-6/N)^2 for every integer N>6','changes_from_R35':changes,'selection_counts':{'strict_improvements':61,'distinct_directions':45,'selected':8},'selected_directions':dirs,'winner_combinations':[{'quadrants':[1],'members':[26,27],'step':'-1/1200000'}],'ledger_counts':{'rows':193,'unique_inputs':193,'new_evaluations':192,'reuse':1},'matches_R26_input_ids':oldmatch,'reused_source_id':'0134','reused_combination_metadata_scope':'row combinations inherited from R35; current task has eight distinct zero-step direction aliases','wall_seconds':str(wall),'child_cpu_seconds':str(cpu),'primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'all_inputs_independently_reconstructed':True,'compressed_archives_match_raw':True,'geometry_or_search_rerun':False,'secant_segments':192,'secant_disjoint_adjacent_pairs':ndisjoint,'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
