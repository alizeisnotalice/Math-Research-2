"""One fixed raw six-atom audit; independent horizontal interval-union subtraction."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,math,os,platform,sys,time
start=time.monotonic();ip=Path(sys.argv[1]);pp=Path(sys.argv[2]);op=Path(sys.argv[3])
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(ip)=='fd212238abff89584129f3471f6b546148cc47e5adc17eed717df3657ba6f874'
assert sha(pp)=='f22e7216fe6a39e36e991fd94a24b00047320376de257fc258e3ebd553f0c9a3'
raw=json.loads(ip.read_text());pr=json.loads(pp.read_text());pts=[tuple(map(Q,z))for z in raw['positions']];mass=list(map(Q,raw['masses']))
assert raw['S_indices']==[0,1,2,3]and Q(raw['threshold'])==1
w=Q(999,1000);e=1-w;assert pts==[(0,0),(w,0),(0,w),(w,w),(Q(-71,200),0),(0,Q(-171,200))]and mass==[Q(1,4)]*5+[Q(499,1000)]
scale=math.lcm(*[z.denominator for p in pts for z in p],*[z.denominator for z in mass])*(1<<72);assert scale==pr['integer_scale']and raw['bits']==pr['bits']==72
counts={'subset_rows':0,'bbox_closures':0,'checked_catalog_fields':0,'radius_brackets':0,'horizontal_slabs':0,'area_endpoint_fields':0}
def box(ids):return(min(pts[i][0]for i in ids),max(pts[i][0]for i in ids),min(pts[i][1]for i in ids),max(pts[i][1]for i in ids))
def subsetmass(mask):return sum((v for i,v in enumerate(mass)if mask>>i&1),Q(0))
def compare(a,b):
 counts['checked_catalog_fields']+=1
 if isinstance(a,Q):assert a==Q(b),(a,b)
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 else:assert a==b,(a,b)
def closure(b):return sum(1<<i for i,(x,y)in enumerate(pts)if b[0]<=x<=b[1]and b[2]<=y<=b[3])
# Distinct coordinate boxes supply a second, complete closure construction.
intervals=[]
for axis in (0,1):
 z=sorted({p[axis]for p in pts});intervals.append([(a,b)for a in z for b in z if a<=b])
coordinate_closures={closure((a,b,c,d))for a,b in intervals[0]for c,d in intervals[1]}-{0}
assert len(coordinate_closures)==23
rows=[];ri=[];ro=[];fi=[];fo=[];fulli=[];fullo=[];eligible=[];subset_closures=set();adj=[set()for i in range(6)]
for mask in range(1,64):
 ids=[i for i in range(6)if mask>>i&1];b=box(ids);m=subsetmass(mask);span=max(b[1]-b[0],b[3]-b[2]);cl=closure(b);subset_closures.add(cl)
 cap=max(Q(1),span*span);assert m<=cap and subsetmass(cl)<=cap
 actual={'mask':mask,'indices':ids,'mass':m,'bbox':list(b),'span_squared':span*span,'cap':cap,'cap_slack':cap-m,'closure_mask':cl,'eligible':m>span*span}
 if actual['eligible']:
  eligible.append(mask)
  for i in ids:adj[i].update(ids)
  rad2=m*scale*scale/4;k=math.isqrt(rad2.numerator//rad2.denominator);up=k if k*k==rad2 else k+1
  assert k*k<=rad2<(k+1)*(k+1)and rad2<=up*up
  counts['radius_brackets']+=1
  a,b0,c,d=[int(x*scale)for x in b];inner=[b0-k,a+k,d-k,c+k];outer=[b0-up,a+up,d-up,c+up]
  actual.update(radius_lower_units=k,radius_upper_units=up,inner=inner,outer=outer)
  fulli.append(inner);fullo.append(outer)
  if mask&15:fi.append(inner);fo.append(outer)
  else:ri.append(inner);ro.append(outer)
 saved=pr['all_source_subset_catalog'][mask-1];assert set(saved)==set(actual)
 for field,val in actual.items():compare(val,saved[field])
 rows.append(actual);counts['subset_rows']+=1
assert subset_closures==coordinate_closures;counts['bbox_closures']=len(subset_closures)
assert [r['mask']for r in pr['distinct_bbox_closures']]==sorted(subset_closures)
for saved in pr['distinct_bbox_closures']:
 mask=saved['mask'];r=rows[mask-1]
 for k,v in saved.items():compare(r[k],v)
L=max(subsetmass(mask)for mask in eligible);maxima=[mask for mask in eligible if subsetmass(mask)==L]
seen={0};todo=[0]
while todo:
 for i in adj[todo.pop()]-seen:seen.add(i);todo.append(i)
span=max(box(range(6))[1]-box(range(6))[0],box(range(6))[3]-box(range(6))[2]);total=sum(mass);gap=span*span-Q(35,33)*total
assert L==1 and maxima==[15]and len(seen)==6 and len(eligible)==12 and gap==Q(395579,250000)>0
qual=pr['qualification'];assert Q(qual['maximum_eligible_mass'])==L and qual['maximum_eligible_masks']==maxima and qual['eligible_subset_count']==len(eligible)and Q(qual['nonterminal_gap_span2_minus_35over33total'])==gap
for saved,want in [('full_inner_old_involving_rectangles',fi),('full_outer_old_involving_rectangles',fo),('complete_residual_inner_rectangles',ri),('complete_residual_outer_rectangles',ro)]:assert pr[saved]==want
assert len(fulli)==12 and len(fi)==9 and len(ri)==3
# Every old-containing receiver lies in C. Pure-outside full labels cancel
# against exactly the same residual BEFORE directed endpoint enclosures.
half=scale//2;wl=int(w*scale);el=int(e*scale);C=(-half,wl+half,-half,wl+half)
assert all(C[0]<=a<b<=C[1]and C[2]<=c<d<=C[3]for a,b,c,d in fo)
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
def enclosing(positive_i,positive_o,negative_i=(),negative_o=(),region=None,excluded=()):
 low=horizontal(positive_i,negative_o,region,excluded);high=horizontal(positive_o,negative_i,region,excluded);assert 0<=low<=high;return low,high
def check_area(saved,ints):
 for k,z in zip(('lower','upper'),ints):assert Q(saved[k])==Q(z,scale*scale);counts['area_endpoint_fields']+=1
 assert Q(saved['width'])==Q(ints[1]-ints[0],scale*scale)
 return dict(lower=str(Q(ints[0],scale*scale)),upper=str(Q(ints[1],scale*scale)))
whole=enclosing(fi,fo,ri,ro);areas={'whole_deletion':check_area(pr['actual_full_plane_deletion_area'],whole)}
cornerboxes={'SW':(-half,-el,-half,-el),'SE':(wl+el,wl+half,-half,-el),'NW':(-half,-el,wl+el,wl+half),'NE':(wl+el,wl+half,wl+el,wl+half)}
cornerints={}
for name,b in cornerboxes.items():
 assert list(map(Q,pr['corner_regions'][name]['bounds']))==[Q(z,scale)for z in b]
 iv=enclosing(fi,fo,ri,ro,b);cornerints[name]=iv;areas[name]=check_area(pr['corner_regions'][name]['area'],iv)
noncorner=enclosing(fi,fo,ri,ro,excluded=tuple(cornerboxes.values()));areas['noncorner']=check_area(pr['actual_noncorner_deletion_area'],noncorner)
assert whole==tuple(sum(v[j]for v in cornerints.values())+noncorner[j]for j in (0,1))
# Full and residual halos themselves, independently from deletion subtraction.
hfull=enclosing(fulli,fullo);hres=enclosing(ri,ro)
areas['whole_full_halo']={k:str(Q(z,scale*scale))for k,z in zip(('lower','upper'),hfull)}
areas['whole_residual_halo']={k:str(Q(z,scale*scale))for k,z in zip(('lower','upper'),hres)}
assert Q(whole[0],scale*scale)<=Q(whole[1],scale*scale)<Q(7,2)
assert Q(cornerints['SW'][0],scale*scale)>Q(1,8)+Q(3,4)*e*e
for name in ('SE','NW','NE'):assert cornerints[name][0]==cornerints[name][1]and Q(cornerints[name][0],scale*scale)==Q(62001,1000000)
# Exact algebra in three independent quadratic radicals; no decimal decision.
radicands=[Q(1,8),Q(999,4000),Q(749,4000)]
def add(*xs):
 out={}
 for x in xs:
  for k,v in x.items():out[k]=out.get(k,Q(0))+v
 return{k:v for k,v in out.items()if v}
def neg(x):return{k:-v for k,v in x.items()}
def mul(x,y):
 out={}
 for i,a in x.items():
  for j,b in y.items():
   z=a*b
   for k,r in enumerate(radicands):
    if i&j&(1<<k):z*=r
   out[i^j]=out.get(i^j,Q(0))+z
 return{k:v for k,v in out.items()if v}
const=lambda x:{0:Q(x)};a={1:Q(1)};U={2:Q(1)};V={4:Q(1)};d=Q(71,200);D=Q(171,200);v=Q(1,4)
first=add(mul(add(mul(const(2),a),const(-d)),add(a,const(-e))),neg(mul(add(a,const(-d+v)),const(v-e))),mul(add(const(d-e),neg(a)),const(v-e)))
joint=add(mul(add(U,const(-e)),add(mul(const(2),U),const(-D))),neg(mul(add(V,const(-e)),add(mul(const(2),V),const(-D)))))
assert first=={0:Q(182323,500000),1:-Q(171,200)}and joint=={0:Q(1,8),2:-Q(857,1000),4:Q(857,1000)}
radbounds=[(Q(176776,500000),Q(176777,500000)),(Q(499749,1000000),Q(1999,4000)),(Q(432723,1000000),Q(432724,1000000))]
for q,(low,up)in zip(radicands,radbounds):assert 0<low<up and low*low<q<up*up
al,au=radbounds[0];ul,uu=radbounds[1];vl,vu=radbounds[2]
assert e<d-au and d-al<v<al and au<d<vl and vu<ul<uu<Q(1,2)and D-uu>au and D-vu>au and D<2*vl and Q(499,4000)<d*d
bound=Q(244823,500000)-D*au-Q(857,1000)*(uu-vl);assert bound==Q(129915191,10**9)>Q(1,8)
# Fine root-directed evaluation also overlaps the independently rebuilt area.
fine=[]
for q in radicands:
 sq=q*scale*scale;k=math.isqrt(sq.numerator//sq.denominator);j=k if k*k==sq else k+1;fine.append((Q(k,scale),Q(j,scale)))
formula_low=Q(244823,500000)-D*fine[0][1]-Q(857,1000)*(fine[1][1]-fine[2][0]);formula_high=Q(244823,500000)-D*fine[0][0]-Q(857,1000)*(fine[1][0]-fine[2][1])
assert max(formula_low,Q(cornerints['SW'][0],scale*scale))<=min(formula_high,Q(cornerints['SW'][1],scale*scale))
out={'status':'PASS','method':'raw63subset and independent36coordinatebox closure qualification; root squaring; horizontal y-slab merged-interval subtraction; original full/residual halos and complete five-region direct decomposition; exact radical algebra','input_sha256':sha(ip),'producer_result_sha256':sha(pp),'script_sha256':sha(__file__),'counts':counts,'qualification':{'eligible_subsets':len(eligible),'unique_max_mask':15,'maximum_mass':'1','connected':True,'nonterminal_gap':str(gap)},'areas':areas,'exact_SW_formula':'244823/500000-(171/200)sqrt(1/8)-(857/1000)(sqrt(999/4000)-sqrt(749/4000))','strict_SW_rational_lower':str(bound),'formula_directed_interval':[str(formula_low),str(formula_high)],'all_five_region_endpoints_sum_exactly':True,'scope':'One actual non-grid balanced six-source measure; SW deletion strictly above local half-mass (and equal allocated single-corner quota), but full deletion below7/2. No uniform-grid or C2 counterexample, no scan, no old-round rerun.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
op.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','counts':counts,'qualification':out['qualification'],'areas':{k:{q:float(Q(z))for q,z in a.items()}for k,a in areas.items()},'runtime':out['runtime']}))
