#!/usr/bin/env python3
"""One fixed five-source qualification/cell, atomic CDF signs, and thin-band constants."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,sys,os,time,platform,hashlib
start=time.perf_counter();ip=Path(sys.argv[1]);cfg=json.loads(ip.read_text());src=[(tuple(map(F,q['xy'])),F(q['mass']),q['old'])for q in cfg['sources']];N=len(src);assert N==5
oldmask=sum(1<<i for i,q in enumerate(src)if q[2]);assert oldmask==15
mass=lambda mask:sum((q[1]for i,q in enumerate(src)if mask>>i&1),F(0))
def box(mask):
 pts=[q[0]for i,q in enumerate(src)if mask>>i&1];return(min(x for x,y in pts),max(x for x,y in pts),min(y for x,y in pts),max(y for x,y in pts))
def span2(b):return max(b[1]-b[0],b[3]-b[2])**2
closures={};eligible=[];parent=list(range(N))
def find(i):
 while parent[i]!=i:i=parent[i]
 return i
for mask in range(1,1<<N):
 b=box(mask);cl=sum(1<<i for i,q in enumerate(src)if b[0]<=q[0][0]<=b[1]and b[2]<=q[0][1]<=b[3]);closures[cl]=(box(cl),mass(cl))
 if mass(mask)>span2(b):
  eligible.append(mask);ids=[i for i in range(N)if mask>>i&1]
  for i in ids[1:]:parent[find(i)]=find(ids[0])
rows=[]
for mask,(b,m)in sorted(closures.items()):
 cap=max(F(1),span2(b));assert m<=cap
 rows.append({'mask':mask,'bbox':list(map(str,b)),'mass':str(m),'span_squared':str(span2(b)),'cap_slack':str(cap-m),'eligible':m>span2(b)})
L=max(map(mass,eligible));maxmasks=[q for q in eligible if mass(q)==L];assert L==1 and maxmasks==[oldmask] and len({find(i)for i in range(N)})==1
nt=span2(box(31))-F(35,33)*mass(31);assert nt>0
lo,hi=map(F,cfg['receiver_cell']);r=F(cfg['common_radius']);assert -F(1,2)<lo<hi<-F(1,10)
maxdist=lambda q:max(abs(z-t)for z in q for t in [lo,hi])
always=[i for i,q in enumerate(src)if maxdist(q[0])<r];assert always==[0,4]
witness=src[0][1];assert -lo<r<F(1,2) and witness-4*r*r==F(1419,250000)>0
outside=src[4][1];assert outside<4*hi*hi and outside<4*F(1,10)**2
area=(hi-lo)**2;quota=F(1,8)+F(3,4)*F(1,10)**2
assert quota==F(53,400) and area==F(151321,1000000)>quota
# Exact four boundary conventions, including coordinates carrying atoms.
def cdf(u,v,le_x,le_y):
 return sum((m for (x,y),m,isold in src if isold and (x<=u if le_x else x<u)and(y<=v if le_y else y<v)),F(0))
assert cdf(F(1,100),F(1,100),False,False)==F(97,100)
ends=list(map(F,cfg['cdf_endpoints']));err=F(cfg['cdf_uncertainty']);queries=0
for l,u in combinations(ends,2):
 for b,v in combinations(ends,2):
  for closed in [False,True]:
   if closed:terms=[(1,u,v,True,True),(-1,l,v,False,True),(-1,u,b,True,False),(1,l,b,False,False)]
   else:terms=[(1,u,v,False,False),(-1,l,v,True,False),(-1,u,b,False,True),(1,l,b,True,True)]
   actual=sum((m for (x,y),m,isold in src if isold and (l<=x<=u and b<=y<=v if closed else l<x<u and b<y<v)),F(0))
   vals=[(sign,cdf(x,y,ex,ey))for sign,x,y,ex,ey in terms];assert sum(sign*z for sign,z in vals)==actual
   lower=max(F(0),sum(sign*(max(F(0),z-err)if sign>0 else min(F(1),z+err))for sign,z in vals))
   upper=min(F(1),sum(sign*(min(F(1),z+err)if sign>0 else max(F(0),z-err))for sign,z in vals))
   assert lower<=actual<=upper;queries+=1
# No new geometric family is searched: these are arithmetic identities only.
parity=[]
for n in cfg['n_values']:
 e=F(1,n*n);K=n//2;g=F((n-2)**2,4*n**4)+(1-2*e)*K*K*e/2;T=(1-2*e)*K*e/2;B=g+T
 expected=F(1,8)-F(1,8*n*n)-F(1,n**3)+F(5,4*n**4)if n%2 else F(1,8)+F(1,4*n)-F(3,2*n**3)+F(1,n**4)
 assert B==expected
 if n%2:assert B<F(1,8)
 else:assert B>F(1,8)+F(3,4*n**4)
 parity.append({'n':n,'thin_joint_band_upper':str(T),'old_single_plus_thin_upper':str(B),'odd_bound_below_one_eighth':n%2==1})
out={'status':'PASS','actual_source_count':N,'actual_subsets_checked':31,'distinct_bbox_closures':len(rows),'closure_rows':rows,'eligible_subset_count':len(eligible),'maximum_eligible_mass':str(L),'unique_maximum_mask':15,'connected':True,'nonterminal_gap':str(nt),'receiver_cell':{'negative_bounds':[str(lo),str(hi)],'common_open_radius':str(r),'always_captured_source_indices':always,'old_label_mass':str(witness),'old_label_strict_margin':str(witness-4*r*r),'actual_captured_mass':str(mass(17)),'residual_total_mass':str(outside),'minimum_residual_threshold':str(4*hi*hi),'area':str(area),'equal_corner_sufficient_quota':str(quota),'area_excess_above_quota':str(area-quota)},'cdf_open_closed_queries':queries,'cdf_lower_upper_interval_checks':queries,'near_origin_cdf_mass':'97/100','thin_band_parity_rows':parity,'scope':'One fixed legal five-source/general-block fixture, complete31subset/closure qualification, positive-area real deletion lower certificate, atomic CDF endpoint identities and scalar thin-band formulas. No full halo union, general area upper bound, random search or repeated production tile calculation.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='closure_rows'}))
