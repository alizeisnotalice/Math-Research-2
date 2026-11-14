"""Independent bounded R13 mass-domination audit; no large halo recomputation.
Checks saved 204 labels, saved exact areas, diff replay, and independent tiny
geometry through the already separate x-column/interval-merging verifier.
"""
from fractions import Fraction as F
from pathlib import Path
import copy,hashlib,importlib.util,json,platform,re,subprocess,sys,time
start=time.time()
implroot=Path(sys.argv[1]); proofroot=Path(sys.argv[2]); evidence=Path(sys.argv[3]); temp=Path(sys.argv[4]); out=Path(sys.argv[5])
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
kernel=load('audited_dominated',implroot/'computation/periodic_halo_dominated.py')
ind=load('independent_small_area',implroot/'computation/independent_periodic_halo.py')
assert digest(implroot/'computation/periodic_halo_dominated.py')=='a697093e30172da71528dc3668add00ea4083811cc2b69ef3a0fc9942eded5e3'
assert digest(implroot/'computation/periodic_halo_transport.py')=='491edb2b6e4763dad22214d16032139c63cb0d59d8427f2cd1694e86e71bd5e9'
# Replay all saved unified hunks against the exact previously audited source.
old=(implroot/'computation/periodic_halo_transport.py').read_text().splitlines(keepends=True)
new=(implroot/'computation/periodic_halo_dominated.py').read_text()
patch=(evidence/'transport_to_dominated.diff').read_text().splitlines(keepends=True)
rebuilt=[];cursor=0;i=0;hunks=0
while i<len(patch):
 m=re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@',patch[i])
 if not m:i+=1;continue
 startold=int(m[1])-1;assert startold>=cursor
 rebuilt.extend(old[cursor:startold]);cursor=startold;i+=1;hunks+=1;oldcount=newcount=0
 while i<len(patch) and not patch[i].startswith('@@ '):
  row=patch[i];tag=row[:1]
  if tag in [' ','-']:
   assert old[cursor]==row[1:];cursor+=1;oldcount+=1
  if tag in [' ','+']:rebuilt.append(row[1:]);newcount+=1
  assert tag in [' ','+','-'];i+=1
 assert oldcount==int(m[2] or 1) and newcount==int(m[4] or 1)
rebuilt.extend(old[cursor:]);assert ''.join(rebuilt)==new
transport=''.join(old)
assert transport[transport.index('    D = 1 << bits'):transport.index('\n\ndef main():')]==new[new.index('    D = 1 << bits'):new.index('\n\ndef main():')]

basepath=implroot/'computation/inputs/r008_final_0611.json'; base=json.loads(basepath.read_text())
probe=json.loads((evidence/'probe204.json').read_text()); saved=json.loads((evidence/'dominated_r5.json').read_text()); meta=saved['transport_matching']
assert digest(basepath)==saved['base_sha256']=='6c2e5d2d8fba5063bd45ecd4467f8c74cca5c0edbef958f90146a103029d5910'
assert digest(evidence/'probe204.json')==saved['input_sha256']
assert meta['base_input']==base
bx,by,lam=map(F,[base['Px'],base['Py'],base['lam']]);px,py=map(F,[probe['Px'],probe['Py']])
assert meta['copy_factors']==[2,2] and (px,py)==(2*bx,2*by)
expected=[(a,b,i) for a in range(2) for b in range(2) for i in range(len(base['pts']))]
assert list(map(tuple,meta['labels']))==expected and len(expected)==204
basepts=[(F(x)%bx,F(y)%by,F(w)) for x,y,w in base['pts']];assert all(w>0 for x,y,w in basepts)
kappa=F(0);maxdist=F(0);actualmass=F(0);liftrecords=[]
for row,label,record in zip(probe['pts'],expected,meta['lifts']):
 x,y,w=map(F,row); x%=px;y%=py;a,b,i=label;sx,sy,sw=basepts[i];zx,zy=sx+a*bx,sy+b*by
 assert w>=0; kx=((zx-x)/px+F(1,2)).__floor__();ky=((zy-y)/py+F(1,2)).__floor__()
 dx=x+kx*px-zx;dy=y+ky*py-zy;distance=max(abs(dx),abs(dy))
 assert distance<=F(meta['delta'])
 actualmass+=w;kappa=max(kappa,w/sw);maxdist=max(maxdist,distance)
 computed={'period_shift':[kx,ky],'displacement':[str(dx),str(dy)],'base_mass':str(sw),'target_mass':str(w),'mass_ratio':str(w/sw)}
 assert computed==record;liftrecords.append(computed)
basemass=sum(w for x,y,w in basepts);assert actualmass>0
prooftext=(proofroot/'proofs/supercell-mass-domination.md').read_text()
printed=re.search(r'm=\\frac\{(\d+)\}\s*\{(\d+)\}',prooftext)
assert printed and F(int(printed[1]),int(printed[2]))==basemass
for dd in [F(1,1000),F(1,10000)]:
    tail0=basemass/(bx*by)*(1+(bx+2*dd)/10)*(1+(by+2*dd)/10)
    assert f'{tail0.numerator} / {tail0.denominator}' in prooftext
R=F(meta['tail_radius']);delta=F(meta['delta'])
def G(m,x,y,r,d):return m/(x*y)*(1+(x+2*d)/(2*r))*(1+(y+2*d)/(2*r))
bound=kappa*G(basemass,bx,by,R,delta)
assert bound<lam and str(bound)==meta['dominated_tail_bound']==saved['tail_bound']
assert str(lam-bound)==meta['tail_margin']
assert str(kappa)==meta['kappa'] and str(maxdist)==meta['maximum_torus_distance']
assert str(actualmass)==meta['target_mass']==saved['mass'] and actualmass==4*basemass
for side in ['inner','outer']:
 assert F(saved[side+'_ratio'])==lam*F(saved[side+'_area'])/actualmass
primary=json.loads((evidence/'primary_r10.json').read_text()); independent=json.loads((evidence/'independent_r10.json').read_text())
assert saved['input_sha256']==primary['input_sha256']==independent['input_sha256']
for side,which in [('inner','lower'),('outer','upper')]:
 assert F(saved[side+'_area'])==F(primary[side+'_area'])==F(independent['area_'+which]['fraction'])
 assert F(saved[side+'_ratio'])==F(primary[side+'_ratio'])==F(independent['weak_ratio_'+which]['fraction'])

# Different tiny examples; all halo geometry below is tiny and separately rebuilt.
def ce(q):return -((-q).__floor__())
def periodic_mass(data,cx,cy,r):
 xperiod,yperiod=F(data['Px']),F(data['Py']);value=F(0)
 for x,y,w in data['pts']:
  x,y,w=F(x)%xperiod,F(y)%yperiod,F(w)
  nx=max(0,((cx+r-x)/xperiod).__floor__()-ce((cx-r-x)/xperiod)+1)
  ny=max(0,((cy+r-y)/yperiod).__floor__()-ce((cy-r-y)/yperiod)+1)
  value+=w*nx*ny
 return value

def independent_area(data,bits=8,cap=F(1)):
 X,Y=F(data['Px']),F(data['Py']);atoms={}
 for x,y,w in data['pts']:
  x,y,w=F(x)%X,F(y)%Y,F(w)
  if not w:continue
  for i in range(ce((-cap-x)/X),((X+cap-x)/X).__floor__()+1):
   for j in range(ce((-cap-y)/Y),((Y+cap-y)/Y).__floor__()+1):
    key=(x+i*X,y+j*Y);atoms[key]=atoms.get(key,F(0))+w
 pts=[(x,y,w) for (x,y),w in atoms.items()]
 grid=ind.setup_grid(pts,(X,Y),bits,cap,F(data['lam']));inside,outside,counts=ind.enumerate_boxes(grid)
 return tuple(F(ind.slab_union_area(rects)[0],grid['S']**2) for rects in (inside,outside)),len(pts)

def small(name,target,basis,delta=F(0),copies=(1,1),require=False):
 res=kernel.calculate_dominated(target,basis,copies,delta,bits=8,R=F(1),require_mass_conservation=require)
 area,n=independent_area(target);assert area==(F(res['inner_area']),F(res['outer_area']))
 mass=sum(F(z[2]) for z in target['pts']);lev=F(target['lam'])
 assert F(res['mass'])==mass and F(res['inner_ratio'])==lev*area[0]/mass
 kap=max(F(w)/F(basis['pts'][i%len(basis['pts'])][2]) for i,(x,y,w) in enumerate(target['pts']))
 assert F(res['transport_matching']['kappa'])==kap
 nchecks=0
 for cx in [F(0),F(1,16),F(1),F(2),F(17,8)]:
  for cy in [F(0),F(1,8),F(2)]:
   for r in [F(1,64),F(1,8),F(1,4),F(1,2),F(1),F(3,2)]:
    assert periodic_mass(target,cx,cy,r)<=kap*periodic_mass(basis,cx,cy,r+delta);nchecks+=1
 return {'input':target,'base':basis,'delta':str(delta),'kappa':str(kap),'mass':str(mass),'areas':list(map(str,area)),'ratio':res['inner_ratio'],'conserved':res['transport_matching']['mass_conserved'],'lifted_tiny_atoms':n,'closed_square_domination_checks':nchecks},res

def rejected(name,target,basis,delta=F(0),**kwargs):
 try:kernel.calculate_dominated(target,basis,(1,1),delta,bits=8,R=F(1),**kwargs)
 except ValueError as e:return {'case':name,'reason':str(e)}
 raise AssertionError(name+' unexpectedly accepted')

smallbase={'Px':'2','Py':'2','lam':'3/2','pts':[['0','0','1/64'],['1/8','0','3/64']]}
zero_collision={'Px':'2','Py':'2','lam':'3/2','pts':[['0','0','3/64'],['0','0','0']]}
tiny={};tiny['zero_target_nonunit_nonconserved_collision'],res=small('zero-collision',zero_collision,smallbase,F(1,8))
assert F(res['transport_matching']['kappa'])==3 and not res['transport_matching']['mass_conserved']
scaledbase=copy.deepcopy(smallbase);scaledbase['lam']='7/2';scaled=copy.deepcopy(zero_collision);scaled['lam']='7/2'
for row in scaled['pts']:row[2]=str(F(row[2])*F(7,3))
tiny['nonunit_7over3_homogeneity'],scaledres=small('scaled',scaled,scaledbase,F(1,8))
assert (res['inner_area'],res['outer_area'],res['inner_ratio'],res['outer_ratio'])==(scaledres['inner_area'],scaledres['outer_area'],scaledres['inner_ratio'],scaledres['outer_ratio'])
collision=copy.deepcopy(zero_collision);collision['pts'][1][2]='1/64'
tiny['two_positive_collision'],positive=small('positive-collision',collision,smallbase,F(1,8),require=True)
assert positive['transport_matching']['kappa']=='3' and positive['transport_matching']['mass_conserved']
tiny['require_mass_conservation_reject']=rejected('mass-flag',zero_collision,smallbase,F(1,8),require_mass_conservation=True)
zero=copy.deepcopy(zero_collision)
for row in zero['pts']:row[2]='0'
tiny['all_zero_reject']=rejected('all-zero',zero,smallbase,F(1,8))
# Exact selected strict-tail boundary at nonunit lambda.
onebase={'Px':'2','Py':'2','lam':'3/2','pts':[['0','0','1/16']]}
near=copy.deepcopy(onebase);near['pts'][0][2]='383/256'
tiny['below_kappa_boundary'],_=small('near-boundary',near,onebase)
for name,w in [('equal','3/2'),('above','385/256')]:
 target=copy.deepcopy(onebase);target['pts'][0][2]=w;tiny[name+'_kappa_boundary_reject']=rejected(name,target,onebase)
# Before-collision maximum is necessary for this labelled certificate, even
# when replacing it by an aggregate ratio would pass the weaker arithmetic.
unequal={'Px':'2','Py':'2','lam':'3/2','pts':[['0','0','1/4096'],['1/8','0','1/16']]}
large={'Px':'2','Py':'2','lam':'3/2','pts':[['0','0','1/2'],['0','0','0']]}
truekap=F(2048);aggregatekap=F(1,2)/sum(F(z[2]) for z in unequal['pts']);baseG=G(sum(F(z[2]) for z in unequal['pts']),F(2),F(2),F(1),F(1,8))
assert truekap*baseG>F(3,2)>aggregatekap*baseG
tiny['collision_large_label_kappa_reject']=rejected('large-label-kappa',large,unequal,F(1,8))
tiny['collision_large_label_kappa_reject'].update({'correct_kappa':str(truekap),'alternative_aggregate_ratio':str(aggregatekap),'correct_tail':str(truekap*baseG),'alternative_aggregate_tail':str(aggregatekap*baseG)})
# Scalar checker generic positivity regression (recorded; repaired script must reject).
scalar=proofroot/'computation/check_supercell_mass_domination_tail.py'
parameter_regressions={}
for name,pxbad,pybad,lambad in [('negative_period_x','-1','1','1'),('negative_period_y','1','-1','1'),('zero_lambda','1','1','0'),('negative_lambda','1','1','-1')]:
 bad={'Px':pxbad,'Py':pybad,'lam':lambad,'pts':[['0','0','1']]};badpath=temp/(name+'.json');badpath.write_text(json.dumps(bad))
 run=subprocess.run([sys.executable,'-B',str(scalar),'--input',str(badpath),'--radius','1','--delta','0','--kappa','1'],capture_output=True,text=True)
 assert run.returncode!=0 and '"status": "passed"' not in run.stdout
 parameter_regressions[name]={'returncode':run.returncode,'stdout_passed':False,'stderr_tail':run.stderr[-200:]}
# Verify monotonicity on exact selected radii, including the true 204-label contract.
rs=[R,R+F(1,1000),R+1,2*R];tails=[kappa*G(basemass,bx,by,r,delta) for r in rs]
assert all(tails[i]>tails[i+1] for i in range(len(tails)-1)) and tails[0]<lam
result={'scope':'R13 bounded independent review; saved 204-label records and tiny geometry, no large geometry rerun',
 'python':platform.python_version(),'script_sha256':digest(Path(__file__)),'diff_hunks_replayed':hunks,
 'probe204':{'input_sha256':digest(evidence/'probe204.json'),'base_sha256':digest(basepath),'rows':204,'copy_factors':[2,2],'kappa':str(kappa),'delta':str(delta),'actual_mass':str(actualmass),'all_saved_lifts_exactly_match':True,'tail_bound':str(bound),'tail_margin':str(lam-bound),'saved_R5_R10_primary_independent_area_and_ratio_equal':True,'independent_lift_records':liftrecords},
 'tiny':tiny,'scalar_parameter_regressions':parameter_regressions,'proof_mass_and_two_tail_table_entries_match_exact_input':True,'geometry_body_after_grid_setup_unchanged':True,
 'tail_monotone_radii':list(map(str,rs)),'tail_monotone_values':list(map(str,tails)),
 'sources':{str(p):digest(p) for p in [implroot/'computation/periodic_halo_dominated.py',implroot/'computation/periodic_halo_transport.py',implroot/'computation/independent_periodic_halo.py',proofroot/'proofs/supercell-mass-domination.md',scalar,evidence/'dominated_r5.json',evidence/'primary_r10.json',evidence/'independent_r10.json']},
 'elapsed_seconds':time.time()-start}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['tiny','probe204']},indent=2));print(json.dumps({'tiny_cases':len(tiny),'probe_rows':204,'all_positive_tiny_geometries_independently_equal':True},indent=2))
