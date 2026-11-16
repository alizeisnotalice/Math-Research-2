"""R81 complete affine caps and two saved endpoints, no parameter search."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import json,hashlib,math,sys,time,platform,os
start=time.monotonic();prod=Path(sys.argv[1]);outpath=Path(sys.argv[2]);own=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(own/'computation/review_r080_balanced_six_atoms.py')=='41dc4d3076c691ca26381d1e88460a499339838133a081c7f538ad9762d0d51e'
# The following four horizontal helpers are copied verbatim from that frozen source.
counts={'horizontal_slabs':0,'subset_boxes':0,'catalog_fields':0,'root_brackets':0,'endpoint_fields':0,'affine_rows':0}
def clip(r,b):
 if b is None:return r
 z=(max(r[0],b[0]),min(r[1],b[1]),max(r[2],b[2]),min(r[3],b[3]));return z if z[0]<z[1]and z[2]<z[3]else None

def merged(iv):
 out=[]
 for a,b in sorted(iv):
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
  else:out.append((a,b))
 return out

def minus_length(pos,neg):
 ans=0
 for a,b in pos:
  left=a
  for c,d in neg:
   if d<=left:continue
   if c>=b:break
   if c>left:ans+=c-left
   left=max(left,d)
   if left>=b:break
  if left<b:ans+=b-left
 return ans

def horizontal(positive,negative,region=None,excluded=()):
 p=[z for r in positive if (z:=clip(r,region))is not None];n=[z for r in negative if (z:=clip(r,region))is not None]+list(excluded)
 if not p:return 0
 yy=sorted({z for r in p+n for z in r[2:]});area=0
 for a,b in zip(yy,yy[1:]):
  y=a+b;pos=merged((l,r)for l,r,c,d in p if 2*c<y<2*d)
  if not pos:continue
  neg=merged((l,r)for l,r,c,d in n if 2*c<y<2*d)
  area+=(b-a)*minus_length(pos,neg);counts['horizontal_slabs']+=1
 return area

fp=prod/'computation/inputs/r081-compute/family.json';sp=prod/'computation/results/round081-compute/production/summary.json';di=prod/'computation/results/round081-compute/regional-credit-diagnostics.json'
FAMILY='e1c332c16b05f0a19171997a79bea6c08690dd43c7204d844e6f435db6b01402';assert sha(fp)==FAMILY
family=json.loads(fp.read_text());summary=json.loads(sp.read_text());diagnostics=json.loads(di.read_text())
p=[tuple(map(Q,z))for z in family['positions']];aff=[tuple(map(Q,z))for z in family['mass_affine_slope_intercept']];w=Q(999,1000);e=Q(1,1000)
assert len(p)==len(aff)==12 and family['S_indices']==[0,1,2,3]and Q(family['w'])==w and Q(family['e'])==e and Q(family['threshold'])==1
assert p[:4]==[(0,0),(w,0),(w,w),(0,w)]and aff[:4]==[(0,Q(1,4))]*4
turn=lambda z:(w-z[1],z[0]);q=(Q(-71,200),0);z=(0,Q(-171,200))
for i in range(4):
 assert p[4+2*i:6+2*i]==[q,z]
 assert aff[4+2*i:6+2*i]==([(0,Q(1,4)),(1,0)]if i%2==0 else[(-1,Q(1,2)),(0,Q(1,4))])
 q,z=turn(q),turn(z)
def bbox(ids):return(min(p[i][0]for i in ids),max(p[i][0]for i in ids),min(p[i][1]for i in ids),max(p[i][1]for i in ids))
def span2(b):return max(b[1]-b[0],b[3]-b[2])**2
def capture(b,ids):return sum(1<<i for i in ids if b[0]<=p[i][0]<=b[1]and b[2]<=p[i][1]<=b[3])
def catalogue(ids):
 ids=list(ids);masks=set()
 for sub in range(1,1<<len(ids)):
  jj=[ids[k]for k in range(len(ids))if sub>>k&1];masks.add(capture(bbox(jj),ids));counts['subset_boxes']+=1
 return [{'mask':mask,'indices':[i for i in ids if mask>>i&1],'bbox':bbox([i for i in ids if mask>>i&1])}for mask in sorted(masks)]
fullcat=catalogue(range(12));rescat=catalogue(range(4,12));oldcat=catalogue(range(4));assert list(map(len,[fullcat,rescat,oldcat]))==[157,104,9]
xx=sorted({x for x,y in p});yy=sorted({y for x,y in p});coordinate={capture((a,b,c,d),range(12))for a in xx for b in xx if a<=b for c in yy for d in yy if c<=d}-{0}
assert len(xx)==len(yy)==6 and len(coordinate)==157 and coordinate=={r['mask']for r in fullcat}
interval=summary['cap_parameter_interval'];rawlo=None;rawhi=None
assert len(interval['rows'])==157
for row,saved in zip(fullcat,interval['rows']):
 ids=row['indices'];s=sum(aff[i][0]for i in ids);b=sum(aff[i][1]for i in ids);cap=max(Q(1),span2(row['bbox']))
 assert saved['mask']==row['mask']and list(map(Q,saved['bbox']))==list(row['bbox'])and Q(saved['mass_slope'])==s and Q(saved['mass_intercept'])==b and Q(saved['cap'])==cap
 if s>0:
  bound=(cap-b)/s;assert saved['direction']=='t_upper'and Q(saved['bound'])==bound;rawhi=bound if rawhi is None else min(rawhi,bound)
 elif s<0:
  bound=(cap-b)/s;assert saved['direction']=='t_lower'and Q(saved['bound'])==bound;rawlo=bound if rawlo is None else max(rawlo,bound)
 else:assert saved['direction']=='constant'and Q(saved['slack'])==cap-b>=0
 counts['affine_rows']+=1
assert rawlo==Q(interval['raw_cap_lower'])==0 and rawhi==Q(interval['raw_cap_upper'])==Q(1,2)
assert list(map(Q,family['requested_t_interval']))==list(map(Q,interval['feasible_intersection']))==[Q(1,4),Q(499,1000)]and not interval['infeasible_constant_masks']
# Raw cap interval endpoints permit zero atoms; the only halo inputs below
# have strictly positive masses. Zero support points do not create mass.
assert all(s*t+b>=0 for s,b in aff for t in [rawlo,rawhi])
def compare(a,b):
 counts['catalog_fields']+=1
 if isinstance(a,Q):assert a==Q(b)
 elif isinstance(a,(list,tuple)):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 else:assert a==b,(a,b)
def root_catalog(cat,masses,L,saved):
 rr=[]
 assert len(cat)==len(saved)
 for r,savedrow in zip(cat,saved):
  ids=r['indices'];m=sum(masses[i]for i in ids);b=r['bbox'];eligible=m>span2(b)
  item={'mask':r['mask'],'indices':ids,'bbox':list(b),'mass':m,'span_squared':span2(b),'eligible':eligible}
  if eligible:
   q=m*L*L/4;lo=math.isqrt(q.numerator//q.denominator);up=lo if lo*lo==q else lo+1;assert lo*lo<=q<(lo+1)**2 and q<=up*up
   a,b0,c,d=[int(z*L)for z in b];item.update(radius_lower_units=lo,radius_upper_units=up,inner=[b0-lo,a+lo,d-lo,c+lo],outer=[b0-up,a+up,d-up,c+up]);counts['root_brackets']+=1
  assert set(item)==set(savedrow)
  for key,val in item.items():compare(val,savedrow[key])
  rr.append(item)
 return rr
results=[]
for name,t in [('lower',Q(1,4)),('upper',Q(499,1000))]:
 ep=prod/f'computation/inputs/r081-compute/endpoints/{name}.json';rp=prod/f'computation/results/round081-compute/production/{name}-result.json'
 saved=json.loads(rp.read_text());inp=json.loads(ep.read_text());m=[s*t+b for s,b in aff];assert all(z>0 for z in m)
 assert list(map(tuple,map(lambda z:map(Q,z),inp['positions'])))==p and list(map(Q,inp['masses']))==m and Q(inp['t'])==t and inp['family_sha256']==FAMILY
 assert saved['endpoint_input_sha256']==sha(ep)and saved['family_sha256']==FAMILY
 L=math.lcm(*[z.denominator for pt in p for z in pt],*[z.denominator for z in m])*(1<<72);assert L==saved['integer_scale']and saved['bits']==72
 full=root_catalog(fullcat,m,L,saved['full_bbox_catalog']);res=root_catalog(rescat,m,L,saved['residual_bbox_catalog']);old=root_catalog(oldcat,m,L,saved['old_bbox_catalog'])
 assert all(r['mass']<=max(Q(1),r['span_squared'])for r in full)
 eligible=[];adj=[set()for i in p]
 # Independent all-subset maximum/ties, not just closure claims.
 for mask in range(1,4096):
  ids=[i for i in range(12)if mask>>i&1];m0=sum(m[i]for i in ids)
  if m0>span2(bbox(ids)):
   eligible.append((mask,m0))
   for i in ids:adj[i].update(ids)
 maximum=max(v for mask,v in eligible);ties=[mask for mask,v in eligible if v==maximum]
 assert maximum==1 and ties==[15,99,390,1548,2073]==saved['qualification']['all_maximum_eligible_masks']and not saved['qualification']['S_unique_maximum']
 seen={0};stack=[0]
 while stack:
  for i in adj[stack.pop()]-seen:seen.add(i);stack.append(i)
 nt=span2(bbox(range(12)))-Q(35,33)*sum(m);assert len(seen)==12 and nt==Q(saved['qualification']['nonterminal_gap'])==Q(45725491,11000000)>0
 assert sum(m)==Q(saved['qualification']['total_mass'])==3 and saved['qualification']['eligible_bbox_count']==sum(r['eligible']for r in full)==25
 for window in saved['shared_side_unit_windows']:
  b=tuple(map(Q,window['unit_closed_bounds']));ids=[i for i in range(12)if capture(b,range(12))>>i&1];assert ids==window['captured_indices']and sum(m[i]for i in ids)==1
  assert sum(m[i]for i in ids if i<4)==sum(m[i]for i in ids if i>=4)==Q(1,2)
 get=lambda cat,key:[r[key]for r in cat if r['eligible']]
 allfi,allfo=get(full,'inner'),get(full,'outer');fi=get([r for r in full if r['mask']&15],'inner');fo=get([r for r in full if r['mask']&15],'outer');ri,ro=get(res,'inner'),get(res,'outer');ii,io=get(old,'inner'),get(old,'outer');assert ii==io
 for key,want in [('full_old_involving_inner',fi),('full_old_involving_outer',fo),('complete_residual_inner',ri),('complete_residual_outer',ro)]:assert saved[key]==want
 half=L//2;wl=int(w*L);el=int(e*L);C=(-half,wl+half,-half,wl+half);corners={'SW':(-half,-el,-half,-el),'SE':(wl+el,wl+half,-half,-el),'NE':(wl+el,wl+half,wl+el,wl+half),'NW':(-half,-el,wl+el,wl+half)}
 assert list(map(Q,saved['C']))==[Q(z,L)for z in C]and all(C[0]<=a<b<=C[1]and C[2]<=c<d<=C[3]for a,b,c,d in fo)
 def enc(ai,ao,bi=(),bo=(),reg=None,excluded=()):
  low=horizontal(ai,bo,reg,excluded);high=horizontal(ao,bi,reg,excluded);assert 0<=low<=high;return low,high
 def endpoints(s,z):
  for k,q in zip(['lower','upper'],z):assert Q(s[k])==Q(q,L*L);counts['endpoint_fields']+=1
  return {k:str(Q(q,L*L))for k,q in zip(['lower','upper'],z)}
 def intersects(xs,ys):return[z for x in xs for y in ys if (z:=clip(x,y))is not None]
 whole=enc(fi,fo,ri,ro);allareas={'whole_deletion':endpoints(saved['whole_plane_actual_deletion'],whole)};ints={};diag=next(c for c in diagnostics['cases']if c['case_id']==name)
 for region,box in list(corners.items())+[('noncorner',C)]:
  excluded=tuple(corners.values())if region=='noncorner'else();sr=saved['five_regions'][region];assert list(map(Q,sr['bounds']))==[Q(z,L)for z in box]
  D=enc(fi,fo,ri,ro,box,excluded);I=enc(ii,io,reg=box,excluded=excluded);ret=enc(ii,io,ri,ro,box,excluded);G=enc(fi,fo,ri+ii,ro+io,box,excluded)
  # Compute T=I intersect residual directly by rectangle intersections,
  # independently of producer's reversal I-retained.
  T=enc(intersects(ii,ri),intersects(ii,ro),reg=box,excluded=excluded)
  assert I[0]==I[1]and D==(ret[0]+G[0],ret[1]+G[1])and T==(I[0]-ret[1],I[1]-ret[0])
  fields={'D':endpoints(sr['actual_deletion'],D),'I':endpoints(sr['internal_halo'],I),'T':endpoints(sr['old_residual_overlap'],T),'old_retained':endpoints(sr['old_retained_deletion'],ret),'G':endpoints(sr['net_mixed_deletion'],G)}
  dr=diag['regions'][region];HF=enc(allfi,allfo,reg=box,excluded=excluded);R=enc(ri,ro,reg=box,excluded=excluded)
  for vals,sv in [(HF,dr['full_halo_area_interval']),(R,dr['complete_residual_area_interval'])]:
   assert [Q(q,L*L)for q in vals]==list(map(Q,sv));counts['endpoint_fields']+=2
  fields.update(full_halo=list(map(lambda q:str(Q(q,L*L)),HF)),complete_residual=list(map(lambda q:str(Q(q,L*L)),R)))
  size=(box[1]-box[0])*(box[3]-box[2])-sum((b-a)*(d-c)for a,b,c,d in excluded)
  assert Q(dr['region_area'])==Q(size,L*L)
  assert list(map(Q,dr['uncovered_full_area_interval']))==[Q(size-HF[1],L*L),Q(size-HF[0],L*L)]
  assert list(map(Q,dr['actual_deletion_deficit_interval']))==[Q(size-D[1],L*L),Q(size-D[0],L*L)]
  assert sr['D_equals_retained_old_plus_net_mixed_at_both_directed_endpoints']is True
  allareas[region]=fields;ints[region]=D
 corner=sum(v[0]for k,v in ints.items()if k!='noncorner'),sum(v[1]for k,v in ints.items()if k!='noncorner')
 allareas['corner_sum']=endpoints(saved['four_corner_actual_deletion_sum'],corner)
 assert whole==(corner[0]+ints['noncorner'][0],corner[1]+ints['noncorner'][1])and Q(corner[1],L*L)<Q(1,2)+3*e*e and Q(whole[1],L*L)<Q(7,2)
 results.append({'case_id':name,'t':str(t),'input_sha256':sha(ep),'producer_result_sha256':sha(rp),'all_maximum_eligible_masks':ties,'eligible_all_subset_count':len(eligible),'eligible_bbox_count':25,'nonterminal_gap':str(nt),'areas':allareas,'five_region_endpoint_identity':True})
out={'status':'PASS','method':'independent exhaustive subset bbox closures and441coordinateboxes; complete affine capacity intersection; exact root squaring; frozen R80 horizontal interval-union algorithm copied with hash binding; complete endpoint5region I/T/G/deletion/full/residual reconstruction','family_sha256':FAMILY,'production_summary_sha256':sha(sp),'regional_diagnostics_sha256':sha(di),'raw_all_cap_t_interval':[str(rawlo),str(rawhi)],'requested_positive_t_interval':['1/4','499/1000'],'zero_mass_boundary_scope':'At raw t0,.5 the measure is nonnegative with zero support weights; no halo computation performed there. Two evaluated endpoints have all12positive atoms.','cases':results,'counts':counts,'reused_horizontal_source_sha256':'41dc4d3076c691ca26381d1e88460a499339838133a081c7f538ad9762d0d51e','script_sha256':sha(__file__),'scope':'One fixed affine family capacity certificate and only two endpoint halo replays for oldS mask15. All5 maximum ties preserved, not a unique-max claim; no middle-t search or uniform area B&B. Analytic whole-t theorem is a separate proof.','runtime':{'pid':os.getpid(),'python':platform.python_version(),'seconds':time.monotonic()-start}}
outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','counts':counts,'cases':[{'id':c['case_id'],'ties':c['all_maximum_eligible_masks'],'corner_sum':c['areas']['corner_sum']}for c in results],'runtime':out['runtime']}))
