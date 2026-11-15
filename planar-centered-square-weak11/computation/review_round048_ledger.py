"""Independent Fraction reconstruction of R43 selection, all inputs and certificate ledger.
No periodic geometry/kernel evaluation is performed by this review.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
from itertools import combinations,product
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
 ap=argparse.ArgumentParser();ap.add_argument('--topic',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();t=a.topic;raw=a.raw;r=t/'computation/results/round048';assert os.path.ismount('/Volumes/data')
 bp=t/'computation/inputs/r008_final_0611.json';op=t/'computation/inputs/r043_final_0195.json';ip=t/'computation/inputs/r048_final_0010.json';mp=t/'computation/inputs/r014_final_0935.json';fp=t/'computation/inputs/r024_common_family.json';vp=t/'computation/results/round017/search_ledger.json';base,target,final,model,family=map(read,(bp,op,ip,mp,fp));source=read(vp)[0];n=len(base['pts']);orig=weights(model);mass=sum(orig);P=[F(base[k])for k in ('Px','Py')];px,py,lam=map(F,(target['Px'],target['Py'],target['lam']))
 assert sha(bp)=='6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910'and sha(op)=='97c70896d07591e65e690be51b08a02b439f070440d3748aad12d0b2ae140dff'and sha(ip)=='b4e2f0f4d8132741747a4569e727431889dcb521840dd053b9a5f87527275e68'and sha(mp)=='498422cc89eecb00436d38ebf84df0028b5c5521db1d762e904d6226292bdd65'and sha(fp)=='cc43a36cd1fa0eeb486fe251fb6b47eca539abcae9ca402c596f14feea9f2322'
 assert source['manifest_sha256']=='9eb82bebf5c15bb91eabd31d8451c96d4d8fe77e1d41b7d14b2d33cde762cac4'and n==51 and len(orig)==len(target['pts'])==len(final['pts'])==204 and mass==sum(weights(target))==sum(weights(final))==4*sum(weights(base))and lam==1 and(px,py)==tuple(2*v for v in P)
 binaryraw=raw.parent.parent/'round017/binary-mass/manifest.json';binarychecked=False
 if binaryraw.exists():
  assert sha(binaryraw)==source['manifest_sha256']and read(binaryraw)['selected_vectors']==source['selected_vectors'];binarychecked=True
 for corner in family['corners']:assert min(weights(corner))>0 and sum(weights(corner))==mass
 vectors=[]
 for v in source['selected_vectors']:
  al=v['aliases'][0];ids=[al['quadrant']*n+j for j in al['members']];anchor=al['quadrant']*n+al['anchor'];G=sum(orig[j]for j in ids);cap=min(G,orig[anchor]);vv=[F(0)]*204
  for j in ids:vv[j]+=al['sign']*F(al['eps'])*cap*orig[j]/G
  vv[anchor]-=al['sign']*F(al['eps'])*cap
  assert sum(vv)==0 and vv==list(map(F,v['vector']));vectors.append(vv)
 assert len(vectors)==8
 def vals(c):return [w+sum(k*v[j]for k,v in zip(c,vectors))for j,w in enumerate(orig)]
 center=[F(-1,8)]*2+[F(1)]*2+[F(42,5)]*4;assert vals(center)==weights(target)
 assert {tuple(vals([aa,aa,F(1),F(1),bb,bb,bb,bb]))for aa,bb in product((F(-1),F(1)),(F(0),F(10)))}=={tuple(weights(c))for c in family['corners']}
 upper=[F(1001,1000)*w for w in weights(base)]*4
 prior=raw.parent.parent/'round043/cross-group-refine';pm=read(prior/'manifest.json');pr=read(prior/'results.json');pz={z['id']:z for z in pr};cache={task['input_sha256']:(task,pz[task['id']],sha(prior/'results.json'))for task in pm['tasks']};tables=[];allhash=set();newhash=set();stages=[];reuses=[];maxk=F(0);minw=None;files=[bp,op,ip,mp,fp,vp,prior/'manifest.json',prior/'results.json'];drivers=['scan_mass_geometry_r48.py','scan_mass_geometry_r48_extension.py','scan_mass_geometry_r48_final_extension.py'];driverhash=['bfd7936a59fc01a356fbc53e40ed736af9b35b55c24ffdc8902ca747c95c3254','425f6ae62eaaf7d0a4d06699c65c98b55027565e0098ed5dadf9d96f3fd06af8','f4e947bd778a9f0ad4ec7786d1abca0f1ca52839a0580cd5029fb091277db7b5'];masscases={'baseline':norm(target)}
 def checkrow(task,row):
  nonlocal maxk,minw
  data=task['input'];ww=weights(data);assert min(ww)>0 and sum(ww)==mass and all(w<=u for w,u in zip(ww,upper))and all(0<=F(x)<px and 0<=F(y)<py for x,y,w in data['pts'])
  dd=delta(base,data);k=max(w/F(base['pts'][i%n][2])for i,w in enumerate(ww));R=F(row['tail_radius']);bound=k*sum(weights(base))/(P[0]*P[1])*(1+(P[0]+2*dd)/(2*R))*(1+(P[1]+2*dd)/(2*R))
  assert 'error'not in row and row['input_sha256']==task['input_sha256']==digest(data)and F(row['delta'])==dd and F(row['kappa'])==k<=F(1001,1000)and F(row['tail_bound'])==bound<1 and R>=5 and F(row['inner_ratio'])<=F(row['outer_ratio'])and row['elapsed_seconds']>0
  if 'delta'in task:assert dd==F(task['delta'])
  if 'box_counts'in row:
   cc=row['box_counts'];assert cc['boxes']==cc['zero_mass']+cc['span_rejected']+cc['root_checked']
  if not row.get('reused',False):assert R==5
  maxk=max(maxk,k);minw=min(ww)if minw is None else min(minw,min(ww))
 for ix,(ir,jr,newcount,reusecount)in enumerate(((range(-4,5),range(-4,5),80,1),(range(0,9),range(-8,1),56,25),(range(4,13),range(-12,-3),56,25))):
  name=f'batch{ix}';d=raw/name;ar=r/name;m=read(d/'manifest.json');rows=read(d/'results.json');su=read(d/'summary.json');done=[json.loads(x)for x in(d/'completed.jsonl').read_text().splitlines()];expected=[]
  for i,j in product(ir,jr):
   coef=[F(i,200)]*2+[F(0)]*2+[F(j,20)]*4;ww=vals([c+v for c,v in zip(center,coef)]);assert min(ww)>0 and all(w<=u for w,u in zip(ww,upper))and sum(ww)==mass;data=norm(target)
   for point,w in zip(data['pts'],ww):point[2]=str(w)
   expected.append({'id':f'{len(expected):04d}','input_sha256':digest(data),'input':data,'coefficients':[list(map(str,coef))]})
  assert m['tasks']==expected and m['schema']=='mass-response-joint-v1'and m['batch']==ix and m['base_sha256']==sha(bp)and m['target_sha256']==sha(op)and m['bits']==40 and m['initial_radius']==5 and m['workers']==10 and m['rejected']==[]and F(m['delta'])==delta(base,target)
  assert m['scan_sha256']==sha(t/'computation/scan_mass_response.py')=='372e2fbf28d53a71372cc4af51e97181bcd86ccee6e05045422140afe2f72908'and m['dominated_kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')=='a697093e30172da71528dc3668add00ea4083811cc2b69ef3a0fc9942eded5e3'
  se=m['selection'];assert se['center']==['-1/8','42/5']and se['base_mass_model_sha256']==sha(mp)and se['source_family_sha256']==sha(fp)and se['source_binary_sha256']==source['manifest_sha256']and se['driver_sha256']==sha(t/'computation'/drivers[ix])==driverhash[ix]and se['family_corners_exactly_reconstructed']is True and se['positive_and_kappa_cap_rejected']==[]and F(se['kappa_cap'])==F(1001,1000)and se['no_fixed_geometry_envelope_inference']is True
  assert len(rows)==len(done)==len(expected)==su['points']==su['completed']==81 and su['new_calls']==newcount and su['reused']==reusecount and su['errors']==su['rejected']==0 and sorted(done,key=lambda z:z['id'])==rows and len({z['id']for z in rows})==81 and digest(m['tasks'])==m['task_list_sha256']==su['manifest_task_list_sha256']and sha(d/'results.json')==su['results_sha256']
  byid={z['id']:z for z in rows};count=0
  for task in expected:
   row=byid[task['id']];checkrow(task,row);assert row.get('reused',False)==(task['input_sha256']in cache)
   if row.get('reused',False):
    oldtask,oldrow,sourcehash=cache[task['input_sha256']];assert oldtask['input']==task['input']and row=={**oldrow,'id':task['id'],'reused':True,'reused_from_id':oldtask['id'],'reused_results_sha256':sourcehash};reuses.append({'batch':name,'id':task['id'],'source_id':oldtask['id'],'source_results_sha256':sourcehash});count+=1
   else:assert task['input_sha256']not in newhash;newhash.add(task['input_sha256'])
   allhash.add(task['input_sha256'])
  assert count==reusecount
  ranked=sorted(rows,key=lambda z:F(z['inner_ratio']),reverse=True)
  for row,top in zip(ranked[:8],su['top']):assert top=={**row,'coefficients':expected[int(row['id'])]['coefficients']}
  assert ranked[0]['id']=='0072';tables.append((m,rows));cache={task['input_sha256']:(task,byid[task['id']],sha(d/'results.json'))for task in expected}
  if ix==0:masscases['first_grid_winner']=expected[72]['input']
  if ix==2:masscases['final_grid_winner']=expected[72]['input']
  for key in ('manifest','results'):assert gzip.decompress((ar/(key+'.json.gz')).read_bytes())==(d/(key+'.json')).read_bytes()
  for key in ('summary.json','completed.jsonl'):assert sha(ar/key)==sha(d/key)
  stages.append({'batch':name,'rows':81,'new':newcount,'reuse':reusecount,'winner':'0072','wall_seconds':su['wall_seconds'],'child_user_cpu_seconds':su['children_user_cpu_seconds'],'child_system_cpu_seconds':su['children_system_cpu_seconds']});files.extend([t/'computation'/drivers[ix],*[d/k for k in ('manifest.json','results.json','summary.json','completed.jsonl')],*sorted(ar.iterdir())])
 assert len(allhash)==193
 # Cross cache contains all three mass stages in order, retaining inherited reuse fields.
 cache={}
 for ix,(mm,rr)in enumerate(tables):
  bz={z['id']:z for z in rr}
  for task in mm['tasks']:cache[task['input_sha256']]=(task,bz[task['id']],sha(raw/f'batch{ix}/results.json'))
 dirs=[('baseline',[],[],F(0)),('q2_28_minus',[2],[28,29],-F(1,500000)),('q2_28_plus',[2],[28,29],F(1,500000)),('all_26_minus',[0,1,2,3],[26,27],-F(1,500000)),('all_26_plus',[0,1,2,3],[26,27],F(1,500000))];expected=[]
 for mn in ('baseline','first_grid_winner','final_grid_winner'):
  for gn,qs,ids,shift in dirs:
   data=norm(masscases[mn])
   for q in qs:
    for j in ids:
     k=q*n+j;data['pts'][k][0]=str((F(data['pts'][k][0])+shift*P[0])%px)
   expected.append({'id':f'{len(expected):04d}','input_sha256':digest(data),'input':data,'delta':str(delta(base,data)),'aliases':[{'mass':mn,'geometry':gn,'shift':str(shift)}]})
 d=raw/'cross';ar=r/'cross';m=read(d/'manifest.json');rows=read(d/'results.json');su=read(d/'summary.json');done=[json.loads(x)for x in(d/'completed.jsonl').read_text().splitlines()];assert m['tasks']==expected and len(rows)==len(done)==su['points']==su['completed']==15 and su['new_calls']==12 and su['reused']==3 and su['errors']==0 and sorted(done,key=lambda z:z['id'])==rows
 assert m['schema']=='r48-mass-geometry-cross-v1'and m['source_batch_results_sha256']==[sha(raw/f'batch{i}/results.json')for i in range(3)]and m['base_sha256']==sha(bp)and m['target_sha256']==sha(op)and m['kernel_sha256']==sha(t/'computation/periodic_halo_dominated.py')and m['scan_sha256']==sha(t/'computation/scan_mass_geometry_r48_cross.py')=='0a3c87790cd74ecae277b6b8575eb96008801a56855e361abd23ae62aaa7d681'and m['bits']==40 and m['initial_R']==5 and m['workers']==10 and digest(m['tasks'])==m['tasks_sha256']and sha(d/'results.json')==su['results_sha256']
 byid={z['id']:z for z in rows};bycase={}
 for task in expected:
  row=byid[task['id']];checkrow(task,row);al=task['aliases'][0];bycase[al['mass'],al['geometry']]=row;assert row.get('reused',False)==(task['input_sha256']in cache)
  if row.get('reused',False):
   oldtask,oldrow,sourcehash=cache[task['input_sha256']];assert oldtask['input']==task['input']and row=={**oldrow,'id':task['id'],'reused':True,'source_results_sha256':sourcehash};reuses.append({'batch':'cross','id':task['id'],'source_id':oldtask['id'],'source_results_sha256':sourcehash})
  else:assert task['input_sha256']not in newhash;newhash.add(task['input_sha256'])
  allhash.add(task['input_sha256'])
 ranked=sorted(rows,key=lambda z:F(z['inner_ratio']),reverse=True)
 for row,top in zip(ranked[:5],su['top']):assert top=={**row,'aliases':expected[int(row['id'])]['aliases']}
 assert ranked[0]['id']=='0010'and norm(expected[10]['input'])==norm(final)==norm(masscases['final_grid_winner'])and sha(d/'candidate_0010.json')==sha(ip)
 assert len(allhash)==205 and len(newhash)==204 and len(reuses)==54
 for key in ('manifest','results'):assert gzip.decompress((ar/(key+'.json.gz')).read_bytes())==(d/(key+'.json')).read_bytes()
 for key in ('summary.json','completed.jsonl'):assert sha(ar/key)==sha(d/key)
 stages.append({'batch':'cross','rows':15,'new':12,'reuse':3,'winner':'0010','wall_seconds':su['wall_seconds'],'child_user_cpu_seconds':su['child_user_cpu_seconds'],'child_system_cpu_seconds':None});files.extend([t/'computation/scan_mass_geometry_r48_cross.py',*[d/k for k in ('manifest.json','results.json','summary.json','completed.jsonl')],*sorted(ar.iterdir())])
 mixed=[]
 for gn,qs,ids,step in dirs[1:]:
  fg,fb,og,ob=[bycase[key]for key in [('final_grid_winner',gn),('final_grid_winner','baseline'),('baseline',gn),('baseline','baseline')]];low=F(fg['inner_ratio'])-F(fb['outer_ratio'])-F(og['outer_ratio'])+F(ob['inner_ratio']);high=F(fg['outer_ratio'])-F(fb['inner_ratio'])-F(og['inner_ratio'])+F(ob['outer_ratio']);assert low<=high and (high<0 or low>0);mixed.append({'geometry':gn,'lower':str(low),'upper':str(high),'sign':'negative'if high<0 else 'positive'})
 assert [x['sign']for x in mixed]==['negative','positive','positive','negative']
 for mn in masscases:
  for gn,qs,ids,step in dirs[1:]:assert F(bycase[mn,gn]['outer_ratio'])<F(bycase[mn,'baseline']['inner_ratio'])
 ea=read(r/'exact_audit.json');fa=read(r/'family_audit.json');assert ea['records']==258 and ea['unique_inputs']==205 and ea['new_calls']==204 and ea['reused_records']==54 and ea['mixed_differences']==mixed and ea['all_12_geometry_changes_strictly_decrease_at_three_tested_masses']is True and ea['script_sha256']==sha(t/'computation/audit_mass_geometry_r48.py')
 for entry,stage in zip(ea['batches'],stages):
  sd=raw/stage['batch'];assert entry['name']==stage['batch']and(entry['points'],entry['new_calls'],entry['reused'],entry['errors'])==(stage['rows'],stage['new'],stage['reuse'],0)and entry['manifest_sha256']==sha(sd/'manifest.json')and entry['results_sha256']==sha(sd/'results.json')and entry['summary_sha256']==sha(sd/'summary.json')
 assert fa['grid_unique_inputs']==193 and fa['cross_aliases']==15 and fa['all_masses_positive']is True and fa['all_masses_exactly_conserved']is True and fa['grid_coordinates_equal_r43']is True and fa['cross_geometries_reconstructed']is True and F(fa['minimum_grid_mass'])==minw==F(1,1250)and F(fa['max_kappa'])==maxk==F(50043,50000)and fa['script_sha256']==sha(t/'computation/audit_family_geometry_r48.py')
 for entry in read(r/'source-manifest.json')['files']:
  q=Path(entry['source']);ar=r/entry['archive'];assert sha(q)==entry['raw_sha256']and sha(ar)==entry['gzip_sha256']and q.stat().st_size==entry['raw_bytes']and ar.stat().st_size==entry['gzip_bytes']and gzip.decompress(ar.read_bytes())==q.read_bytes()
 p,ind=map(read,(r/'primary_old_R10.json',r/'independent_old_R10.json'));lo,up=map(F,(p['inner_ratio'],p['outer_ratio']));assert p['input_sha256']==ind['input_sha256']==sha(ip)and p['code_sha256']==sha(t/'computation/periodic_halo.py')=='b788671299ee33c14e4c6d872e36bfab06cbaa55dd6cd6a1c9bc6c26190a0f33'and ind['script_sha256']==sha(t/'computation/independent_periodic_halo.py')=='d1e05a0627d699b54472b1b3e74ad500035135853298d33f519f78d3a17903e8'
 assert lo==F(ind['weak_ratio_lower']['fraction'])==F(ranked[0]['inner_ratio'])and up==F(ind['weak_ratio_upper']['fraction'])==F(ranked[0]['outer_ratio'])and F(p['inner_area'])==F(ind['area_lower']['fraction'])==lo*mass and F(p['outer_area'])==F(ind['area_upper']['fraction'])==up*mass
 assert [F(p[k])for k in ('period_x','period_y','level','mass')]==[px,py,lam,mass]and list(map(F,ind['periods']))==[px,py]and F(ind['lambda'])==lam and F(ind['base_mass'])==mass and p['bits']==ind['bits']==40 and F(p['tail_radius'])==F(ind['radius_cap'])==10
 assert p['counts']['boxes']==ind['counts']['all_span_admissible_boxes']and p['window_atoms_distinct']==ind['lifted_atoms']and p['inner_rectangles_unique']==ind['counts']['unique_inner_rectangles']and p['outer_rectangles_unique']==ind['counts']['unique_outer_rectangles']
 tail=mass/(px*py)*(1+px/20)*(1+py/20);q=math.ceil(F(10)/min(px,py))+1;assert tail==F(p['tail_bound'])==F(ind['tail_bound_at_cap']['fraction'])<1 and q==3
 for name in ('primary_old_R10.json','independent_old_R10.json'):assert sha(raw/'cross'/name)==sha(r/name)
 prev=t/'computation/results/round043/primary_old_R10.json';gain=lo-F(read(prev)['outer_ratio']);assert gain>0
 changes=[i for i,(x,y)in enumerate(zip(final['pts'],target['pts']))if F(x[2])!=F(y[2])];assert changes==[1,28,29,32,33,52,79,80,103,130,131,134,135,154,181,182]and all(x[:2]==y[:2]for x,y in zip(final['pts'],target['pts']))and vals([F(-13,200)]*2+[F(1)]*2+[F(39,5)]*4)==weights(final)
 proof=t/'proofs/round048-lower.md';report=t/'rounds/048.md';fr=[F(int(x),int(y))for x,y in re.findall(r'\\frac\{\s*(\d+)\s*\}\s*\{\s*(\d+)\s*\}',proof.read_text())];assert fr==[tail,lo,up,gain]and sha(proof)=='faec25876cadb2753360da45098c3e3f5940ba72de1399fc84b964d88a8350ba'and sha(report)=='fb6044acda90a7bd109c855a5061f2cc6722005503864719f26c864968f97988'
 c=read(r/'integration_checks.json');assert c['input_sha256']==sha(ip)and c['primary_sha256']==sha(r/'primary_old_R10.json')and c['independent_sha256']==sha(r/'independent_old_R10.json')and F(c['lower'])==lo and F(c['upper'])==up and F(c['gain_over_r43_upper'])==gain and F(c['tail_bound'])==tail and c['q']==3 and c['changed_mass_label_indices_zero_based']==changes
 for st in stages:
  cc=c['counts'][st['batch']];assert(cc['points'],cc['new_calls'],cc['reused'])==(st['rows'],st['new'],st['reuse'])and cc['wall_seconds']==st['wall_seconds']
 files.extend([prev,proof,report,t/'proofs/periodic-complete-halo.md',t/'computation/scan_mass_response.py',t/'computation/scan_supercell_mass.py',t/'computation/scan_quadrant_groups.py',t/'computation/scan_supercell_2x2.py',t/'computation/periodic_halo_dominated.py',t/'computation/periodic_halo.py',t/'computation/independent_periodic_halo.py',t/'computation/audit_mass_geometry_r48.py',t/'computation/audit_family_geometry_r48.py',t/'computation/build_round048_ledger.py',*[x for x in r.iterdir()if x.is_file()]])
 out={'status':'passed','scope':'Complete four-batch exact reconstruction, eight mass transfers and family corners, all reuse chains, mixed intervals and frozen two-kernel certificates; no geometry replay','input_sha256':sha(ip),'canonical_input_sha256':digest(final),'mass_parameters':{'a':'-13/200','b':'39/5'},'positions_unchanged_from_R43':True,'changed_mass_indices_from_R43':changes,'eight_zero_sum_mass_vectors_independently_reconstructed':True,'original_R17_binary_manifest_identity_checked':binarychecked,'original_R17_binary_manifest_sha256':source['manifest_sha256'],'family_four_corners_reconstructed':True,'minimum_mass':str(minw),'maximum_kappa':str(maxk),'stages':stages,'ledger_counts':{'rows':258,'unique_inputs':205,'new_evaluations':204,'reuse':54},'reuses':reuses,'mixed_differences':mixed,'all_12_tested_geometry_changes_strictly_decrease':True,'lower':str(lo),'upper':str(up),'lower_minus_R43_full_upper':str(gain),'tail_bound_R10':str(tail),'q':q,'D':str(2**40),'epsilon':str(F(1,2**41)),'finite_L1_factor':'(1-6/N)^2 for every integer N>6','primary_recorded_seconds':p['elapsed_seconds'],'independent_recorded_seconds':ind['total_seconds'],'geometry_or_search_rerun':False,'source_vector_scope':'Uses archived R17 selected-vector ledger; reconstructs every vector from mass group/anchor/epsilon formula, not a fresh R15 search ranking','fixed_geometry_family_scope':'Mass corners are a source identity check; R24 fixed-R14-geometry envelope does not bound the new R43 geometry','files_sha256':{str(x):sha(x)for x in files},'audit_script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('files_sha256','reuses')},indent=2))
if __name__=='__main__':main()
