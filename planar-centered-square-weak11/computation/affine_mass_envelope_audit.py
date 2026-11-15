"""R24 independent tiny full-box envelope and B&B audit.
Direct per-point mass sums, Fraction endpoints, and a simple slab merge; no
prefix table or original rectangle union is used for independent bounds.
No current 204-label geometry or search is executed.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import copy,hashlib,importlib.util,json,math,subprocess,sys,time
start=time.process_time();wall_start=time.monotonic();topic=Path(sys.argv[1]);theory=Path(sys.argv[2]);temp=Path(sys.argv[3]);out=Path(sys.argv[4]);comp=topic/'computation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sys.modules[name]=m;sp.loader.exec_module(m);return m
old=load('periodic_halo',comp/'periodic_halo.py');new=load('review_envelope',comp/'mass_family_envelope.py')
def ce(q):return -((-q).__floor__())
def area(rects):
 if not rects:return F(0)
 xs=sorted({z for r in rects for z in r[:2]});total=F(0)
 for l,h in zip(xs,xs[1:]):
  intervals=sorted((b,t) for x0,x1,b,t in rects if x0<=l<h<=x1)
  length=F(0)
  if intervals:
   a,b=intervals[0]
   for c,d in intervals[1:]:
    if c>b:length+=b-a;a,b=c,d
    else:b=max(b,d)
   length+=b-a
  total+=(h-l)*length
 return total

def point(data,a,b):
 a0,a1=map(F,data['a_bounds']);b0,b1=map(F,data['b_bounds']);x=(a-a0)/(a1-a0);y=(b-b0)/(b1-b0);corners=data['corners']
 result={k:corners[0][k] for k in ['Px','Py','lam']};result['pts']=[]
 for rows in zip(*(c['pts'] for c in corners)):
  w=list(F(r[2]) for r in rows);result['pts'].append([rows[0][0],rows[0][1],str(w[0]+x*(w[1]-w[0])+y*(w[2]-w[0]))])
 return result

def direct(data,cell,bits,R,inner=False):
 pts=point(data,F(cell[0]),F(cell[2]));P,Q,lam=map(F,[pts['Px'],pts['Py'],pts['lam']]);D=1<<bits
 cellcorners=[point(data,a,b) for a,b in product([F(cell[0]),F(cell[1])],[F(cell[2]),F(cell[3])])]
 # List all labelled periodic copies directly. Collisions stay as separate
 # summands; rectangle masses are summed from points, not a prefix table.
 copies=[]
 for i,(x,y,w) in enumerate(pts['pts']):
  x,y=F(x)%P,F(y)%Q
  for j in range(ce((-R-x)/P),((P+R-x)/P).__floor__()+1):
   for k in range(ce((-R-y)/Q),((Q+R-y)/Q).__floor__()+1):copies.append((x+j*P,y+k*Q,i))
 xs=sorted({x for x,y,i in copies});ys=sorted({y for x,y,i in copies});rects=set();boxes=0;cell_empty=0
 for a,b in product(xs,xs):
  if a>b or b-a>2*R:continue
  for c,d in product(ys,ys):
   if c>d or d-c>2*R:continue
   inds=[i for x,y,i in copies if a<=x<=b and c<=y<=d]
   if not inds:continue
   boxes+=1
   masses=[sum(F(p['pts'][i][2]) for i in inds) for p in cellcorners];m=max(masses)
   span=max(b-a,d-c)
   if not m or span>=2*R or lam*span*span>=m:cell_empty+=1;continue
   root_squared=m*D*D/(4*lam);k=math.isqrt(root_squared.numerator//root_squared.denominator)
   assert F(k*k)<=root_squared<F((k+1)*(k+1))
   radius=min(R,F(max(0,k-1) if inner else k+1,D))
   x0,x1=max(F(0),b-radius),min(P,a+radius);y0,y1=max(F(0),d-radius),min(Q,c+radius)
   if x0<x1 and y0<y1:rects.add((x0,x1,y0,y1))
 mass=sum(F(z[2]) for z in pts['pts']);return area(rects),mass,{'all_boxes':boxes,'cell_empty_exact':cell_empty,'rectangles':len(rects)}

def family(kind='critical',scale=F(1),lam=F(1),domain=None):
 base={'Px':'10','Py':'10','lam':str(lam),'pts':[['0','0','1/2'],['1','0','1/2'],['5','5','1']]};cs=[]
 for a,b in [(-1,-1),(1,-1),(-1,1),(1,1)]:
  d=copy.deepcopy(base)
  w=[F(1,2)+F(a,8)+F(b,16),F(1,2)+F(a,8)-F(b,16),(F(2) if kind=='different_base_mass' else F(1))-F(a,4)]
  for i,value in enumerate(w):d['pts'][i][2]=str(scale*value)
  if kind=='collision':d['pts'][1][:2]=['0','0']
  if kind=='wrap':d['pts'][0][0]='-1/4'
  cs.append(d)
 return {'base':base,'copy_factors':[1,1],'delta':'1' if kind=='collision' else '1/4' if kind=='wrap' else '0','a_bounds':(domain[0] if domain else ['-1','1']),'b_bounds':(domain[1] if domain else ['-1','1']),'corners':cs}

cases=[('critical',family(),F(2),[2,5]),('collision',family('collision'),F(2),[2,5]),('wrap',family('wrap'),F(2),[2,5]),('nonunit_lambda',family(scale=F(3,2),lam=F(3,2)),F(2),[2,5]),('actual_mass_not_base',family('different_base_mass'),F(2),[2,5]),('cap_nondyadic_R',family(scale=F(1,2)),F(3,4),[0,3]),('nonunit_parameter_domain',family(domain=([F(-2),F(3)],[F(1),F(4)])),F(2),[2,5])]
records=[];point_checks=0
for name,data,R,bitslist in cases:
 a0,a1=map(F,data['a_bounds']);b0,b1=map(F,data['b_bounds']);am,bm=(a0+a1)/2,(b0+b1)/2
 cells=[(a0,a1,b0,b1),(a0,am,b0,bm),(am,a1,bm,b1),(am,am,b0,b1),(am,am,bm,bm),(a0,a0,b0,b0),(a0+(a1-a0)*F(1,4),a0+(a1-a0)*F(2,3),b0+(b1-b0)*F(3,8),b0+(b1-b0)*F(7,8))]
 for bits in bitslist:
  f=new.MassFamilyEnvelope(data,bits,R);assert f.meta['all_corner_labels_strictly_positive']
  expected_mass=sum(F(z[2]) for z in data['corners'][0]['pts']);assert f.mass==expected_mass
  for cell in cells:
   measured=f.upper(*cell);ind,mass,stats=direct(data,cell,bits,R)
   assert F(measured['upper_area'])==ind and F(measured['upper_ratio'])==ind*F(data['base']['lam'])/mass
   if cell[0]==cell[1] and cell[2]==cell[3]:
    p=point(data,F(cell[0]),F(cell[2]));old_point=old.calculate(p,bits=bits,tail=R)
    assert F(old_point['outer_area'])==ind and F(old_point['outer_ratio'])==F(measured['upper_ratio']);point_checks+=1
   records.append({'family':name,'bits':bits,'R':str(R),'cell':list(map(str,cell)),'upper_area':str(ind),'upper_ratio':measured['upper_ratio'],'direct_stats':stats})
# The exact root .5 at bits2 is deliberately rounded outward by one full grid.
assert F(3,4)-F(1,2)==F(1,4)
critical_result=next(r for r in records if r['family']=='critical' and r['bits']==2 and r['cell']==['0','0','0','0'])
assert critical_result['upper_area']=='17/4'
# Fixed bits mesh-only floor: a dummy-parameter constant family has the same U
# on every cell. True halo area is2 (joint near-point label is exactly inactive).
constant=family();constant['corners']=[copy.deepcopy(constant['base']) for _ in range(4)]
fixed=[]
for bits in [2,5,8]:
 f=new.MassFamilyEnvelope(constant,bits,F(2));u=f.upper(-1,1,-1,1);v=f.upper(0,0,0,0)
 assert u['upper_area']==v['upper_area'] and F(u['upper_area'])>2
 fixed.append({'bits':bits,'all_mesh_cells_and_point_upper_area':u['upper_area'],'true_area':'2'})
assert all(F(fixed[i]['all_mesh_cells_and_point_upper_area'])>F(fixed[i+1]['all_mesh_cells_and_point_upper_area']) for i in range(2))
# Generic clipped-rectangle increment bound, including old empty/inverted boxes.
band_checks=0;P,Q=F(3,2),F(5,4)
for box in [(F(0),F(0),F(0),F(0)),(F(0),F(1),F(0),F(1)),(F(-2),F(-1),F(1,2),F(2)),(F(1),F(3),F(-1),F(1))]:
 def rect_area(radius):
  a,b,c,d=box;return max(F(0),min(P,a+radius)-max(F(0),b-radius))*max(F(0),min(Q,c+radius)-max(F(0),d-radius))
 for radius,eps in product([F(0),F(1,10),F(1,2),F(1)],[F(1,100),F(1,5)]):
  difference=rect_area(radius+eps)-rect_area(radius)
  assert 0<=difference<=2*eps*(P+Q)+4*eps*eps;band_checks+=1

# Real B&B/pilot tiny run, then independent leaf-cover and bound replay.
positive=family();inp=temp/'positive-family.json';inp.write_text(json.dumps(positive))
def run_driver(script,args,filename,expect_success=True):
 output=temp/filename;cmd=[sys.executable,'-B',str(comp/script),str(inp),'--output',str(output),*args]
 run=subprocess.run(cmd,capture_output=True,text=True)
 if not expect_success:
  assert run.returncode!=0;return {'returncode':run.returncode,'stderr_tail':run.stderr[-250:]}
 assert run.returncode==0,run.stderr
 return json.loads(output.read_text())
bb=run_driver('branch_mass_family.py',['--bits','5','--radius','2','--max-leaves','4','--seconds','30','--gap','1/10000000000'],'bb.json')
leaves=bb['leaves'];assert len(leaves)<=4 and len(leaves)==bb['divisions']+1
cells=[tuple(map(F,r['cell'])) for r in leaves];assert sum((ah-al)*(bh-bl) for al,ah,bl,bh in cells)==4
for i,(al,ah,bl,bh) in enumerate(cells):
 assert F(-1)<=al<ah<=1 and -1<=bl<bh<=1
 for cl,ch,dl,dh in cells[i+1:]:assert min(ah,ch)<=max(al,cl) or min(bh,dh)<=max(bl,dl)
for leaf,cell in zip(leaves,cells):
 u,m,_=direct(positive,cell,5,F(2));assert F(leaf['upper'])==u/m
 if leaf['status']=='pruned':assert F(leaf['upper'])<F(bb['lower'])
for sample in bb['samples']:
 a,b=map(F,sample['point']);low,m,_=direct(positive,(a,a,b,b),5,F(2),inner=True);assert low/m==F(sample['lower'])
assert F(bb['lower'])==max(F(p['lower']) for p in bb['samples'])
assert F(bb['upper'])==max([F(bb['lower'])]+[F(l['upper']) for l in leaves])
if bb['status']=='target_gap':assert F(bb['gap'])<=F(bb['target_gap'])
pilot=run_driver('pilot_mass_family.py',['--bits','3','--radius','2','--workers','1','--reference-lower','0'],'pilot.json')
last=None
for grid in pilot['grids']:
 side=math.isqrt(grid['partitions']);assert side*side==grid['partitions']==len(grid['cells'])
 cell_list=[tuple(map(F,c['cell'])) for c in grid['cells']];assert sum((a1-a0)*(b1-b0) for a0,a1,b0,b1 in cell_list)==4
 for cell,record in zip(cell_list,grid['cells']):
  u,m,_=direct(positive,cell,3,F(2));assert F(record['upper'])==u/m
 value=max(F(c['upper']) for c in grid['cells']);assert value==F(grid['upper']) and (last is None or value<=last);last=value
# Zero-corner mode is valid for upper bounds but excludes the strict rate budget.
zero=copy.deepcopy(positive)
for j,(u,v) in enumerate([(F(0),F(0)),(F(1),F(0)),(F(0),F(1)),(F(1),F(1))]):
 weights=[u/8,F(1,2)+v/16,F(3,2)-u/8-v/16]
 for row,w in zip(zero['corners'][j]['pts'],weights):row[2]=str(w)
zero_kernel=new.MassFamilyEnvelope(zero,3,F(2));assert not zero_kernel.meta['all_corner_labels_strictly_positive']
zero_input=temp/'zero-family.json';zero_input.write_text(json.dumps(zero))
run=subprocess.run([sys.executable,'-B',str(comp/'branch_mass_family.py'),str(zero_input),'--output',str(temp/'zero-rejected.json'),'--bits','3','--radius','2','--max-leaves','1'],capture_output=True,text=True)
assert run.returncode!=0 and 'strictly positive' in run.stderr
zero_point=point(zero,F(-1),F(-1));zero_old=old.calculate(zero_point,bits=3,tail=F(2))
assert zero_kernel.upper(-1,-1,-1,-1)['upper_area']==zero_old['outer_area']
run=subprocess.run([sys.executable,'-B',str(comp/'branch_mass_family.py'),str(zero_input),'--output',str(temp/'zero-allowed.json'),'--bits','3','--radius','2','--max-leaves','1','--allow-zero-corners'],capture_output=True,text=True)
assert run.returncode==0,run.stderr
zero_mode=json.loads((temp/'zero-allowed.json').read_text());assert not zero_mode['family']['all_corner_labels_strictly_positive']
rejects={}
def reject(name,data):
 try:new.MassFamilyEnvelope(data,3,F(2))
 except ValueError as error:rejects[name]=str(error);return
 raise AssertionError(name+' accepted')
bad=copy.deepcopy(positive);bad['corners'][0]['pts'][0][2]='-1';reject('negative_corner',bad)
bad=copy.deepcopy(positive);bad['corners'][0]['pts'][0][2]='1/3';reject('nonaffine_corner',bad)
bad=copy.deepcopy(positive)
for j in [1,3]:bad['corners'][j]['pts'][2][2]=str(F(bad['corners'][j]['pts'][2][2])+F(1,10))
reject('affine_but_nonconstant_total',bad)
bad=copy.deepcopy(positive);bad['corners'][0]['pts'][0][0]='1/100';reject('varying_geometry',bad)
bad=family(scale=F(10));reject('uniform_tail_fails',bad)
# Helper exactly reproduces the toy two-box area; no project halo is inferred.
helper=theory/'computation/check_affine_mass_enclosure.py';run=subprocess.run([sys.executable,'-B',str(helper)],capture_output=True,text=True,check=True);helper_result=json.loads(run.stdout);assert helper_result['status']=='passed' and helper_result['union_area']=='9'
result={'scope':'independent tiny full-box and driver audit only; no real 204-label evaluation',
 'python':sys.version,'main_process_cpu_seconds':time.process_time()-start,'elapsed_wall_seconds':time.monotonic()-wall_start,'script_sha256':sha(Path(__file__)),
 'independent_envelope_comparisons':len(records),'point_cell_old_outer_equalities':point_checks,'cases':records,
 'regression_initial_point_mismatch':{'before_fix_new_area':'9/2','old_and_repaired_area':'17/4','cause':'globally possible box exactly inactive at the point; corrected by per-cell exact span prune'},
 'fixed_bits_mesh_floor':fixed,'clipped_rectangle_increment_checks':band_checks,
 'tiny_branch_and_bound':bb,'tiny_pilot':pilot,'zero_corner_default_driver_rejected':True,'zero_corner_explicit_mode_runs_without_strict_rate_claim':True,'zero_corner_point_old_outer_equal':True,'invalid_family_rejections':rejects,
 'helper_actual_output':helper_result,'sources':{str(p):sha(p) for p in [comp/'mass_family_envelope.py',comp/'branch_mass_family.py',comp/'pilot_mass_family.py',comp/'periodic_halo.py',theory/'proofs/affine-mass-rectangle-enclosure.md',helper]}}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','tiny_branch_and_bound','tiny_pilot']},indent=2))
