"""Independent Fraction reconstruction of R35 two grids, all inputs, reuse, and secants.
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
 ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round035';assert os.path.ismount('/Volumes/data')
 bp=t/'computation/inputs/r008_final_0611.json';op=t/'computation/inputs/r029_final_0298.json';ip=t/'computation/inputs/r035_final_0134.json';base,old,final=map(read,(bp,op,ip));n=len(base['pts']);P=[F(base[k])for k in ('Px','Py')];px,py,lam=map(F,(final['Px'],final['Py'],final['lam']));mass=sum(weights(final));kappa=max(w/F(base['pts'][i%n][2])for i,w in enumerate(weights(old)))
 assert sha(bp)=='6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910' and sha(op)=='ca9f8b0ea3d78036151ed3e5c37ffd558c395fc2de66c273f561f4ef86630f59'and sha(ip)=='5021c54704af45b3ef8289c67a245926473434640372a934b2f15c4485f0c026'
 assert n==51 and len(final['pts'])==len(old['pts'])==204 and min(weights(final))>0 and weights(final)==weights(old)and mass==4*sum(weights(base))and(px,py)==tuple(2*v for v in P)and lam==1 and kappa==F(25021,25000)
 assert all(F(old[k])==F(final[k])for k in ('Px','Py','lam'))and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in final['pts'])
 changes={name:[i for i,(u,v)in enumerate(zip(old['pts'],final['pts']))if F(u[axis])!=F(v[axis])]for axis,name in enumerate(('x','y','mass'))};assert changes=={'x':[130,131],'y':[],'mass':[]}
 assert all(F(final['pts'][i][0])-F(old['pts'][i][0])==F(2958614645661,250000000000000000)==P[0]*F(9,2500000)for i in changes['x'])
 prior=raw.parent/'round026/quadrant-groups';prior29=raw.parent/'round029/quadrant-joint';pm=read(prior/'manifest.json');pr=read(prior/'results.json');s29=read(prior29/'summary.json');z29=read(prior29/'results.json');pt={z['id']:z for z in pm['tasks']};pz={z['id']:z for z in pr}
 assert sha(prior/'manifest.json')=='4445cfc88d7ba633bd8b3d495a1aa90a119d542a751b23fe917931e85953f691'and sha(prior/'results.json')=='0e7633d75509ef95a353a05901d45b8a3605c4199e725d28e93d6974fd7759fb'
 improved=sorted((z for z in pr if F(z['inner_ratio'])>F(pz['0000']['outer_ratio'])),key=lambda z:F(z['inner_ratio']),reverse=True);keys={(pt[z['id']]['aliases'][0]['quadrant'],pt[z['id']]['aliases'][0]['group']['axis'],tuple(pt[z['id']]['aliases'][0]['group']['members']),F(pt[z['id']]['aliases'][0]['shift'])>0)for z in improved};assert len(improved)==61 and len(keys)==45
 directions=[{'quadrants':[q]}for q in range(4)]+[{'quadrants':[0,1,2,3]}];allhash=set();newhash=set();stages=[];reuses=[];secants=[];files=[bp,op,ip,prior/'manifest.json',prior/'results.json',prior29/'summary.json',prior29/'results.json'];cache={};ctask={}
 for ix,(name,limit,den,points,newcount,reusecount,script,codehash)in enumerate((('axis-wall-probe',40,100000,401,400,1,'scan_axis_wall_probe_r35.py','248bb79bef06bfbc9451283dfe6717207bf544c3054f6e2c9d64e44d16873b52'),('axis-wall-refine',25,2500000,251,240,11,'scan_axis_wall_refine_r35.py','07473f72902a7eeaa1a02bcb7ba33fd7e442333bf96485f4e5c5cfe17bd54426'))):
  d=raw/name;ar=r/name;m=read(d/'manifest.json');rows=read(d/'results.json');su=read(d/'summary.json');done=[json.loads(x)for x in(d/'completed.jsonl').read_text().splitlines()];expected={}
  for direction in directions:
   for j in range(-limit,limit+1):
    data=norm(old)
    for q in direction['quadrants']:
     for label in (28,29):
      k=q*n+label;data['pts'][k][0]=str((F(data['pts'][k][0])+F(j,den)*P[0])%px)
    h=digest(data)
    if h not in expected:expected[h]={'id':'','input_sha256':h,'input':data,'delta':str(delta(base,data)),'combinations':[]}
    expected[h]['combinations'].append({'quadrants':direction['quadrants'],'step':str(F(j,den))})
  h=digest(old);baseline=expected.pop(h);expected={h:baseline,**expected}
  for i,z in enumerate(expected.values()):z['id']=f'{i:04d}'
  assert m['tasks']==list(expected.values())and len(expected)==len(rows)==len(done)==su['points']==su['completed']==points and su['new_calls']==newcount and su['reused']==reusecount and su['errors']==0
  assert m['schema']=='r35-single-axis-probe-v1'and m['directions']==directions and m['strict_improving_count']==61 and m['distinct_direction_count']==45 and m['base_sha256']==sha(bp)and m['target_sha256']==sha(op)and m['source_manifest_sha256']==sha(prior/'manifest.json')and m['source_results_sha256']==sha(prior/'results.json')and m['source_r29_summary_sha256']==sha(prior29/'summary.json')
  assert m['kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')=='a697093e30172da71528dc3668add00ea4083811cc2b69ef3a0fc9942eded5e3'and m['script_sha256']==sha(t/'computation'/script)==codehash and m['bits']==40 and m['radius']==8 and m['workers']==12 and digest(m['tasks'])==m['task_list_sha256']and sha(d/'results.json')==su['results_sha256']
  assert len({z['id']for z in rows})==len({z['input_sha256']for z in rows})==points and sorted(done,key=lambda z:z['id'])==rows
  byid={z['id']:z for z in rows};count=0
  for task in m['tasks']:
   row=byid[task['id']];data=task['input'];dd=delta(base,data);bound=kappa*sum(weights(base))/(P[0]*P[1])*(1+(P[0]+2*dd)/16)*(1+(P[1]+2*dd)/16)
   assert weights(data)==weights(old)and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in data['pts'])and 'error'not in row and row['input_sha256']==task['input_sha256']and F(row['delta'])==dd==F(task['delta'])and F(row['kappa'])==kappa and F(row['tail_radius'])==8 and F(row['tail_bound'])==bound<1 and F(row['inner_ratio'])<=F(row['outer_ratio'])
   cc=row['box_counts'];assert cc['boxes']==cc['zero_mass']+cc['span_rejected']+cc['root_checked']and 0<=cc['inner_empty']<=cc['root_checked']and 0<=cc['outer_empty']<=cc['root_checked']
   shouldreuse=(task['id']=='0000')if ix==0 else task['input_sha256']in cache
   assert row.get('reused',False)==shouldreuse
   if shouldreuse:
    if ix==0:
     oldrow=s29['top'][0];source=prior29/'results.json';assert oldrow['id']=='0298'and {k:v for k,v in oldrow.items()if k!='combinations'}==next(z for z in z29 if z['id']=='0298')and oldrow['input_sha256']==digest(old)
    else:
     oldrow=cache[task['input_sha256']];source=raw/'axis-wall-probe/results.json';assert ctask[task['input_sha256']]['input']==data and ctask[task['input_sha256']]['delta']==task['delta']
    assert row=={**oldrow,'id':task['id'],'reused':True,'source_results_sha256':sha(source)};reuses.append({'batch':name,'id':task['id'],'source_id':oldrow['id'],'source_results_sha256':sha(source)});count+=1
   else:assert task['input_sha256']not in newhash;newhash.add(task['input_sha256'])
   allhash.add(task['input_sha256'])
  assert count==reusecount
  order=sorted(rows,key=lambda z:F(z['inner_ratio']),reverse=True)
  for row,top in zip(order[:8],su['top']):assert top=={**row,'combinations':m['tasks'][int(row['id'])]['combinations']}
  win=order[0]
  assert win['id']==('0000'if ix==0 else '0134')
  if ix==1:assert norm(m['tasks'][134]['input'])==norm(final)and sha(d/'candidate_0134.json')==sha(ip)and m['tasks'][134]['combinations']==[{'quadrants':[2],'step':'9/2500000'}]
  # Independently form all exact slope intervals, including baseline aliases.
  analysis=read(d/'secant_analysis.json');groups={}
  for task in m['tasks']:
   for alias in task['combinations']:groups.setdefault(tuple(alias['quadrants']),[]).append((F(alias['step']),byid[task['id']]))
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
  assert analysis['directions']==expected_analysis and analysis['manifest_sha256']==sha(d/'manifest.json')and analysis['results_sha256']==sha(d/'results.json')and analysis['script_sha256']==sha(t/'computation/analyze_axis_probe_r35.py')
  for fname in ('manifest','results'):assert gzip.decompress((ar/(fname+'.json.gz')).read_bytes())==(d/(fname+'.json')).read_bytes()
  for fname in ('summary.json','completed.jsonl','secant_analysis.json'):assert sha(ar/fname)==sha(d/fname)
  stages.append({'batch':name,'rows':points,'new':newcount,'reuse':reusecount,'aliases':sum(len(x['combinations'])for x in m['tasks']),'winner':win['id'],'winner_lower':win['inner_ratio'],'wall_seconds':su['wall_seconds'],'child_cpu_seconds':str(Decimal(str(su['child_user_cpu_seconds']))+Decimal(str(su['child_system_cpu_seconds'])))})
  secants.append({'batch':name,'sample_points_per_direction':2*limit+1,'segments':5*2*limit,'disjoint_adjacent_pairs':ndisjoint})
  if ix==0:cache={z['input_sha256']:z for z in rows};ctask={z['input_sha256']:z for z in m['tasks']}
  files.extend([t/'computation'/script,*[d/x for x in ('manifest.json','results.json','summary.json','completed.jsonl','secant_analysis.json')],*sorted(ar.iterdir())])
 assert len(allhash)==641 and len(newhash)==640 and len(reuses)==12
 p,ind=map(read,(r/'primary_old_R10.json',r/'independent_old_R10.json'));lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));assert p['input_sha256']==ind['input_sha256']==sha(ip)
 assert p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
 assert lo==F(ind['weak_ratio_lower']['fraction'])==F(stages[1]['winner_lower'])and up==F(ind['weak_ratio_upper']['fraction'])and lo<=up and F(p['inner_area'])==F(ind['area_lower']['fraction'])==mass*lo and F(p['outer_area'])==F(ind['area_upper']['fraction'])==mass*up
 assert [F(p[k])for k in ('period_x','period_y','level','mass')]==[px,py,lam,mass]and list(map(F,ind['periods']))==[px,py]and F(ind['lambda'])==lam and F(ind['base_mass'])==mass and p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==10
 tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])<1 and q==3
 assert p['counts']['boxes']==ind['counts']['all_span_admissible_boxes']and p['window_atoms_distinct']==ind['lifted_atoms']and p['inner_rectangles_unique']==ind['counts']['unique_inner_rectangles']and p['outer_rectangles_unique']==ind['counts']['unique_outer_rectangles']
 prev=t/'computation/results/round029/primary_old_R10.json';gain=lo-F(read(prev)['outer_ratio']);assert gain>0
 proof=t/'proofs/round035-lower.md';report=t/'rounds/035.md';fr=[F(int(x),int(y))for x,y in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]and sha(proof)=='50389429d03645402801e03e6a27f62af7b4c3af7b00991dd7066ac801a6fe3a'and sha(report)=='784e21375ff2afe3fe9826291d303b2c00f82692aa16f3faba25e784b1e7227f'
 c=read(r/'integration_checks.json');assert c['input_sha256']==sha(ip)and c['primary_sha256']==sha(r/'primary_old_R10.json')and c['independent_sha256']==sha(r/'independent_old_R10.json')and F(c['lower'])==lo and F(c['upper'])==up and F(c['gain_over_r29_upper'])==gain and F(c['tail_bound'])==tail and c['q']==3 and c['changed_label_indices_zero_based']==changes['x']
 for stage in stages:
  sc=c['counts'][stage['batch']];su=read(raw/stage['batch']/'summary.json');assert(sc['points'],sc['new_calls'],sc['reused'])==(stage['rows'],stage['new'],stage['reuse'])and sc['wall_seconds']==stage['wall_seconds']and sc['child_cpu_seconds']==su['child_user_cpu_seconds']+su['child_system_cpu_seconds']
  assert all(sha(path)in report.read_text()for path in [raw/stage['batch']/'manifest.json',raw/stage['batch']/'results.json',r/stage['batch']/'manifest.json.gz',r/stage['batch']/'results.json.gz'])
 for name in ('primary_old_R10.json','independent_old_R10.json'):assert sha(raw/'axis-wall-refine'/name)==sha(r/name)
 files.extend([r/'integration_checks.json',t/'computation/build_round035_ledger.py'])
 files.extend([prev,proof,report,t/'proofs/periodic-complete-halo.md',t/'computation/periodic_halo.py',t/'computation/independent_periodic_halo.py',t/'computation/periodic_halo_dominated.py',t/'computation/analyze_axis_probe_r35.py',r/'primary_old_R10.json',r/'independent_old_R10.json'])
 out={'status':'passed','scope':'Both complete ledgers, exact secant enclosures and frozen complete-period certificates; no geometry replay','input_sha256':sha(ip),'canonical_input_sha256':digest(final),'lower':str(lo),'upper':str(up),'lower_minus_R29_full_upper':str(gain),'tail_bound_R10':str(tail),'q':q,'D':str(2**40),'epsilon':str(F(1,2**41)),'finite_L1_factor':'(1-6/N)^2 for every integer N>6','changes_from_R29':changes,'stages':stages,'ledger_counts':{'rows':652,'unique_inputs':641,'new_evaluations':640,'reuse':12},'reuses':reuses,'secant_analysis':secants,'secant_scope':'Disjoint adjacent enclosures certify non-affinity on their combined sampled interval; neither exact wall nor continuous monotonicity is certified','primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'all_inputs_independently_reconstructed':True,'compressed_archives_match_raw':True,'geometry_or_search_rerun':False,'files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='files_sha256'},indent=2))
if __name__=='__main__':main()
