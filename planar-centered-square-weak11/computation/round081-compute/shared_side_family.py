#!/usr/bin/env python3
"""One exact alternating family: cap interval and its two feasible endpoints."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,os,platform,sys,time
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
ROOT=Path(__file__).resolve().parents[2];CORE=ROOT/'computation/round075/certify_three_atom_collapse.py'
assert sha(CORE)=='a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456'
spec=importlib.util.spec_from_file_location('frozen_r075',CORE);area=importlib.util.module_from_spec(spec);spec.loader.exec_module(area)
def closures(points,indices):
 xx=sorted({points[j][0] for j in indices});yy=sorted({points[j][1] for j in indices})
 def intervals(axis,values):
  groups=[sum(1<<j for j in indices if points[j][axis]==x) for x in values];out=[]
  for i in range(len(groups)):
   z=0
   for j in range(i,len(groups)):z|=groups[j];out.append(z)
  return out
 masks={x&y for x in intervals(0,xx) for y in intervals(1,yy)}-{0};out=[]
 for mask in sorted(masks):
  jj=[j for j in indices if mask>>j&1];box=(min(points[j][0] for j in jj),max(points[j][0] for j in jj),min(points[j][1] for j in jj),max(points[j][1] for j in jj));span=max(box[1]-box[0],box[3]-box[2]);out.append({'mask':mask,'indices':jj,'bbox':box,'span_squared':span*span})
 return out,{'axis_coordinate_counts':[len(xx),len(yy)],'candidate_axis_rectangle_count':len(xx)*(len(xx)+1)*len(yy)*(len(yy)+1)//4,'distinct_nonempty_bbox_closures':len(out)}
def cap_interval(catalog,aff,requested):
 lower=upper=None;bad=[];rows=[]
 for row in catalog:
  slope=sum(aff[j][0] for j in row['indices']);intercept=sum(aff[j][1] for j in row['indices']);cap=max(F(1),row['span_squared']);z={'mask':row['mask'],'bbox':list(map(str,row['bbox'])),'mass_slope':str(slope),'mass_intercept':str(intercept),'cap':str(cap)}
  if slope>0:
   bound=(cap-intercept)/slope;upper=bound if upper is None else min(upper,bound);z.update(direction='t_upper',bound=str(bound))
  elif slope<0:
   bound=(cap-intercept)/slope;lower=bound if lower is None else max(lower,bound);z.update(direction='t_lower',bound=str(bound))
  else:
   z.update(direction='constant',slack=str(cap-intercept))
   if intercept>cap:bad.append(row['mask'])
  rows.append(z)
 lo=max(requested[0],lower) if lower is not None else requested[0];hi=min(requested[1],upper) if upper is not None else requested[1]
 return {'raw_cap_lower':None if lower is None else str(lower),'raw_cap_upper':None if upper is None else str(upper),'requested_interval':list(map(str,requested)),'feasible_intersection':list(map(str,[lo,hi])) if not bad and lo<=hi else None,'infeasible_constant_masks':bad,'rows':rows}
def directed(catalog,masses,L):
 out=[]
 for z in catalog:
  mass=sum(masses[j] for j in z['indices']);r={'mask':z['mask'],'indices':z['indices'],'bbox':list(map(str,z['bbox'])),'mass':str(mass),'span_squared':str(z['span_squared']),'eligible':mass>z['span_squared']}
  if r['eligible']:
   k2=mass*L*L;assert k2.denominator==1;k=math.isqrt(k2.numerator//4);up=k if 4*k*k==k2 else k+1;assert 4*k*k<=k2<=4*up*up
   a,b,c,d=[int(v*L) for v in z['bbox']];r.update(radius_lower_units=k,radius_upper_units=up,inner=[b-k,a+k,d-k,c+k],outer=[b-up,a+up,d-up,c+up])
  out.append(r)
 return out
def endpoint(case_id,t,p,masses,w,e,bits,fullcat,rescat,oldcat,outdir,inputdir,familysha):
 started=time.monotonic();assert all(m>0 for m in masses);L=math.lcm(*[z.denominator for q in p for z in q],*[m.denominator for m in masses])*(1<<bits)
 raw={'schema':'r081-alternating-actual-endpoint-v1','threshold':'1','t':str(t),'positions':[list(map(str,q)) for q in p],'masses':list(map(str,masses)),'S_indices':[0,1,2,3],'family_sha256':familysha};ip=inputdir/f'{case_id}.json';assert not ip.exists();ip.write_text(json.dumps(raw,indent=2)+'\n')
 full=directed(fullcat,masses,L);res=directed(rescat,masses,L);old=directed(oldcat,masses,L)
 for r in full:assert F(r['mass'])<=max(F(1),F(r['span_squared']))
 elig=[r for r in full if r['eligible']];maxmass=max(F(r['mass']) for r in elig);ties=[r['mask'] for r in elig if F(r['mass'])==maxmass];assert maxmass==1 and 15 in ties
 parent=list(range(12))
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for r in elig:
  for j in r['indices'][1:]:parent[find(j)]=find(r['indices'][0])
 span=max(max(x for x,y in p)-min(x for x,y in p),max(y for x,y in p)-min(y for x,y in p));nt=span*span-F(35,33)*sum(masses)
 def rr(cat,key):return [r[key] for r in cat if r['eligible']]
 oldfull=[r for r in full if r['mask']&15];fi,fo=rr(oldfull,'inner'),rr(oldfull,'outer');ri,ro=rr(res,'inner'),rr(res,'outer');ii,io=rr(old,'inner'),rr(old,'outer');assert ii==io
 half=L//2;wl=int(w*L);el=int(e*L);C=[-half,wl+half,-half,wl+half];Q={'SW':[-half,-el,-half,-el],'SE':[wl+el,wl+half,-half,-el],'NE':[wl+el,wl+half,wl+el,wl+half],'NW':[-half,-el,wl+el,wl+half]};assert all(C[0]<=a<b<=C[1] and C[2]<=c<d<=C[3] for a,b,c,d in fo)
 def crop(rects,box):
  if box is None:return rects
  out=[]
  for a,b,c,d in rects:
   z=[max(a,box[0]),min(b,box[1]),max(c,box[2]),min(d,box[3])]
   if z[0]<z[1] and z[2]<z[3]:out.append(z)
  return out
 def bound(FI,FO,RI,RO,box=None,exclude=()):
  x,y=crop(FI,box),crop(FO,box);u,v=crop(RI,box)+list(exclude),crop(RO,box)+list(exclude);lo=area.event_area(x,v);hi=area.event_area(y,u);lo2,c1=area.midpoint_area(x,v);hi2,c2=area.midpoint_area(y,u);assert (lo,hi)==(lo2,hi2) and 0<=lo<=hi
  return {'lower':str(F(lo,L*L)),'upper':str(F(hi,L*L)),'lower_decimal':float(F(lo,L*L)),'upper_decimal':float(F(hi,L*L)),'two_algorithms_identical':True,'midpoint_cells':[c1,c2]},(lo,hi)
 quota=F(1,2)+3*e*e;regions={};total,tint=bound(fi,fo,ri,ro);sumlo=sumhi=0
 for name,box in list(Q.items())+[('noncorner',C)]:
  exclusion=list(Q.values()) if name=='noncorner' else [];D,ints=bound(fi,fo,ri,ro,box,exclusion);I,iints=bound(ii,io,[],[],box,exclusion);ret,rtints=bound(ii,io,ri,ro,box,exclusion);G,gints=bound(fi,fo,ri+ii,ro+io,box,exclusion)
  assert ints==(rtints[0]+gints[0],rtints[1]+gints[1]);assert iints[0]==iints[1]
  # I is exact for these four fixed old corners; old-residual overlap
  # endpoints are the reversals of retained-old endpoints.
  overlap={'lower':str(F(iints[0]-rtints[1],L*L)),'upper':str(F(iints[1]-rtints[0],L*L))}
  reg={'bounds':list(map(lambda x:str(F(x,L)),box)),'actual_deletion':D,'internal_halo':I,'old_residual_overlap':overlap,'old_retained_deletion':ret,'net_mixed_deletion':G,'D_equals_retained_old_plus_net_mixed_at_both_directed_endpoints':True}
  if name!='noncorner':
   sumlo+=ints[0];sumhi+=ints[1];reg['geometric_quota_quarter']=str(quota/4);reg['D_minus_quota_interval']=[str(F(D['lower'])-quota/4),str(F(D['upper'])-quota/4)]
  regions[name]=reg
 nc=regions['noncorner']['actual_deletion'];assert tint==(sumlo+int(F(nc['lower'])*L*L),sumhi+int(F(nc['upper'])*L*L))
 corner={'lower':str(F(sumlo,L*L)),'upper':str(F(sumhi,L*L)),'lower_decimal':float(F(sumlo,L*L)),'upper_decimal':float(F(sumhi,L*L))};sidewindows=[];unitbox=(F(-171,200),1-F(171,200),F(0),F(1))
 for i in range(4):
  near=4+2*i;far=5+2*((i-1)%4);assert F(1,2)+masses[near]+masses[far]==1
  captured=[j for j,(x,y) in enumerate(p) if unitbox[0]<=x<=unitbox[1] and unitbox[2]<=y<=unitbox[3]];assert len(captured)==4 and near in captured and far in captured and sum(masses[j] for j in captured)==1
  sidewindows.append({'corner_near_index':near,'previous_corner_far_index':far,'unit_closed_bounds':list(map(str,unitbox)),'captured_indices':captured,'old_edge_mass':'1/2','outside_sum':str(masses[near]+masses[far]),'unit_capture_mass':'1'})
  a,b,c,d=unitbox;unitbox=(w-d,w-c,a,b)
 result={'schema':'r081-shared-side-endpoint-result-v1','pid':os.getpid(),'source_sha256':sha(__file__),'family_sha256':familysha,'endpoint_input_sha256':sha(ip),'case_id':case_id,'t':str(t),'bits':bits,'integer_scale':L,'full_bbox_catalog':full,'residual_bbox_catalog':res,'old_bbox_catalog':old,'qualification':{'all_bbox_caps_pass':True,'maximum_eligible_mass':str(maxmass),'all_maximum_eligible_masks':ties,'S_mask':15,'S_unique_maximum':ties==[15],'eligible_bbox_count':len(elig),'connected_eligible_hyperedges':len({find(j) for j in range(12)})==1,'nonterminal_gap':str(nt),'nonterminal':nt>0,'total_mass':str(sum(masses))},'shared_side_unit_windows':sidewindows,'full_old_involving_inner':fi,'full_old_involving_outer':fo,'complete_residual_inner':ri,'complete_residual_outer':ro,'C':list(map(lambda x:str(F(x,L)),C)),'whole_plane_actual_deletion':total,'five_regions':regions,'four_corner_actual_deletion_sum':corner,'aggregate_corner_quota':str(quota),'aggregate_corner_upper_pass':F(corner['upper'])<=quota,'aggregate_corner_lower_violation':F(corner['lower'])>quota,'aggregate_quota_minus_corner_interval':[str(quota-F(corner['upper'])),str(quota-F(corner['lower']))],'whole_plane_upper_below_3point5':F(total['upper'])<F(7,2),'two_algorithms':'vertical event/interval union and compressed Cartesian midpoint coverage, exact integer endpoints','pureoutside_full_residual_cancelled_exactly':True,'scope':'One feasible endpoint of the frozen shared-side family. All corner quantities use the same outside measure. Per-corner excess/deficit are bookkeeping; only aggregate is tested. No uniform parameter claim.','seconds':time.monotonic()-started}
 op=outdir/f'{case_id}-result.json';assert not op.exists();op.write_text(json.dumps(result,indent=2)+'\n')
 return {'case_id':case_id,'t':str(t),'input':str(ip),'input_sha256':sha(ip),'result':op.name,'result_sha256':sha(op),'qualification':result['qualification'],'whole_plane_actual_deletion':total,'four_corner_actual_deletion_sum':corner,'corner_values':{j:regions[j]['actual_deletion'] for j in Q},'aggregate_quota':str(quota),'aggregate_corner_upper_pass':result['aggregate_corner_upper_pass'],'aggregate_corner_lower_violation':result['aggregate_corner_lower_violation'],'seconds':result['seconds']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--family',required=True);ap.add_argument('--output-dir',required=True);ap.add_argument('--endpoint-input-dir',required=True);a=ap.parse_args();started=time.monotonic();assert Path('/Volumes/data').is_mount();fp=Path(a.family);data=json.loads(fp.read_text());assert data['schema']=='r081-shared-side-alternating-family-v1';p=[tuple(map(F,z)) for z in data['positions']];aff=[tuple(map(F,z)) for z in data['mass_affine_slope_intercept']];assert len(p)==12 and data['S_indices']==[0,1,2,3];w=F(data['w']);e=F(data['e']);assert w==F(999,1000) and e==1-w
 fullcat,fullstats=closures(p,range(12));rescat,resstats=closures(p,range(4,12));oldcat,oldstats=closures(p,range(4));assert fullstats['candidate_axis_rectangle_count']==441
 interval=cap_interval(fullcat,aff,list(map(F,data['requested_t_interval'])));outdir=Path(a.output_dir);inputdir=Path(a.endpoint_input_dir);outdir.mkdir(parents=True,exist_ok=True);inputdir.mkdir(parents=True,exist_ok=True);assert not (outdir/'summary.json').exists()
 summary={'schema':'r081-shared-side-family-summary-v1','pid':os.getpid(),'python':sys.version,'platform':platform.platform(),'source_sha256':sha(__file__),'family_input_sha256':sha(fp),'frozen_union_core_sha256':sha(CORE),'full_menu_statistics':fullstats,'residual_menu_statistics':resstats,'old_menu_statistics':oldstats,'cap_parameter_interval':interval,'endpoints':[],'uniform_parameter_certificate_performed':False,'parameter_search_performed':False}
 if interval['feasible_intersection'] is not None:
  ts=list(map(F,interval['feasible_intersection']))
  for case_id,t in [('lower',ts[0])]+([('upper',ts[1])] if ts[1]!=ts[0] else []):summary['endpoints'].append(endpoint(case_id,t,p,[s*t+b for s,b in aff],w,e,int(data['root_bits']),fullcat,rescat,oldcat,outdir,inputdir,sha(fp)))
 summary['seconds']=time.monotonic()-started;(outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:summary[k] for k in ['pid','full_menu_statistics','residual_menu_statistics','endpoints','uniform_parameter_certificate_performed','seconds']},indent=2));print('cap interval',interval['raw_cap_lower'],interval['raw_cap_upper'],interval['feasible_intersection'])
if __name__=='__main__':main()
