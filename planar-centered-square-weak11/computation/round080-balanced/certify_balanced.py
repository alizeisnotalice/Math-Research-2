#!/usr/bin/env python3
"""One actual six-source measure: complete caps, qualification and exact areas."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,os,platform,sys,time
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
ROOT=Path(__file__).resolve().parents[2];CORE=ROOT/'computation/round075/certify_three_atom_collapse.py'
assert sha(CORE)=='a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456'
spec=importlib.util.spec_from_file_location('frozen_r075',CORE);area=importlib.util.module_from_spec(spec);spec.loader.exec_module(area)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);args=ap.parse_args();start=time.monotonic();ip=Path(args.input);data=json.loads(ip.read_text());p=[tuple(map(F,q)) for q in data['positions']];m=list(map(F,data['masses']));w=F(data['w']);e=F(data['e']);bits=int(data['bits'])
 assert p==[(F(0),F(0)),(w,F(0)),(F(0),w),(w,w),(F(-71,200),F(0)),(F(0),F(-171,200))] and m==[F(1,4)]*5+[F(499,1000)] and w==F(999,1000) and e==1-w and data['S_indices']==[0,1,2,3] and F(data['threshold'])==1
 L=math.lcm(*[z.denominator for q in p for z in q],*[z.denominator for z in m])*(1<<bits);rows=[];closures={};parent=list(range(6))
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 def join(i,j):parent[find(j)]=find(i)
 for mask in range(1,64):
  jj=[j for j in range(6) if mask>>j&1];mass=sum(m[j] for j in jj);box=(min(p[j][0] for j in jj),max(p[j][0] for j in jj),min(p[j][1] for j in jj),max(p[j][1] for j in jj));span=max(box[1]-box[0],box[3]-box[2]);cap=max(F(1),span*span);eligible=mass>span*span
  closure=sum(1<<j for j,(x,y) in enumerate(p) if box[0]<=x<=box[1] and box[2]<=y<=box[3]);cmass=sum(m[j] for j in range(6) if closure>>j&1);assert cmass<=cap and mass<=cap
  closures[closure]={'mask':closure,'bbox':list(map(str,box)),'mass':str(cmass),'span_squared':str(span*span),'cap':str(cap),'cap_slack':str(cap-cmass),'eligible':cmass>span*span}
  row={'mask':mask,'indices':jj,'mass':str(mass),'bbox':list(map(str,box)),'span_squared':str(span*span),'cap':str(cap),'cap_slack':str(cap-mass),'closure_mask':closure,'eligible':eligible}
  if eligible:
   for j in jj[1:]:join(jj[0],j)
   z=mass*L*L;assert z.denominator==1;k=math.isqrt(z.numerator//4);up=k if 4*k*k==z else k+1;assert 4*k*k<=z<=4*up*up
   a,b,c,d=[int(v*L) for v in box];row.update(radius_lower_units=k,radius_upper_units=up,inner=[b-k,a+k,d-k,c+k],outer=[b-up,a+up,d-up,c+up])
  rows.append(row)
 eligible=[r for r in rows if r['eligible']];maxmass=max(F(r['mass']) for r in eligible);maximum=[r['mask'] for r in eligible if F(r['mass'])==maxmass];components={find(j) for j in range(6)}
 totalmass=sum(m);spanall=max(max(x for x,y in p)-min(x for x,y in p),max(y for x,y in p)-min(y for x,y in p));nonterminalgap=spanall*spanall-F(35,33)*totalmass
 residual=[r for r in rows if not r['mask']&15];oldfull=[r for r in rows if r['mask']&15];assert len(rows)==63 and len(residual)==3 and len(oldfull)==60
 # Full pure-outside labels cancel EXACTLY against the identical residual
 # halo before directed bounds, avoiding artificial residual boundary strips.
 def rr(source,key):return [r[key] for r in source if r['eligible']]
 fi,fo,ri,ro=rr(oldfull,'inner'),rr(oldfull,'outer'),rr(residual,'inner'),rr(residual,'outer')
 C=[int(-F(1,2)*L),int((w+F(1,2))*L),int(-F(1,2)*L),int((w+F(1,2))*L)]
 assert all(C[0]<=a<b<=C[1] and C[2]<=c<d<=C[3] for a,b,c,d in fo)
 def clip(rectangles,box):
  if box is None:return rectangles
  out=[]
  for a,b,c,d in rectangles:
   z=[max(a,box[0]),min(b,box[1]),max(c,box[2]),min(d,box[3])]
   if z[0]<z[1] and z[2]<z[3]:out.append(z)
  return out
 def enc(box=None,excluded=()):
  fulli,fullo=clip(fi,box),clip(fo,box);resi=clip(ri,box)+list(excluded);reso=clip(ro,box)+list(excluded)
  lo=area.event_area(fulli,reso);hi=area.event_area(fullo,resi);lo2,nlo=area.midpoint_area(fulli,reso);hi2,nhi=area.midpoint_area(fullo,resi);assert (lo,hi)==(lo2,hi2) and 0<=lo<=hi
  return {'lower':str(F(lo,L*L)),'upper':str(F(hi,L*L)),'lower_decimal':float(F(lo,L*L)),'upper_decimal':float(F(hi,L*L)),'width':str(F(hi-lo,L*L)),'two_algorithms_identical':True,'midpoint_cells':[nlo,nhi]},(lo,hi)
 wl=int(w*L);el=int(e*L);half=L//2;cornerboxes={'SW':[-half,-el,-half,-el],'SE':[wl+el,wl+half,-half,-el],'NW':[-half,-el,wl+el,wl+half],'NE':[wl+el,wl+half,wl+el,wl+half]}
 whole,whint=enc();corners={};sumlo=sumhi=0
 for name,box in cornerboxes.items():
  value,ints=enc(box);corners[name]={'bounds':list(map(lambda z:str(F(z,L)),box)),'area':value};sumlo+=ints[0];sumhi+=ints[1]
 noncorner,ncint=enc(excluded=list(cornerboxes.values()));assert (whint[0],whint[1])==(sumlo+ncint[0],sumhi+ncint[1])
 sw=corners['SW']['area'];halfP=F(1,8);main=F(1,8)+F(3,4)*e*e
 out={'schema':'r080-balanced-six-atom-certificate-v1','pid':os.getpid(),'python':sys.version,'platform':platform.platform(),'input_sha256':sha(ip),'source_sha256':sha(__file__),'frozen_area_core_sha256':sha(CORE),'bits':bits,'integer_scale':L,'source_count':6,'full_subset_count':63,'old_involving_subset_count':60,'residual_subset_count':3,'all_source_subset_catalog':rows,'distinct_bbox_closures':[closures[k] for k in sorted(closures)],'distinct_bbox_closure_count':len(closures),'all_subset_and_bbox_caps_pass':True,'qualification':{'maximum_eligible_mass':str(maxmass),'maximum_eligible_masks':maximum,'S_mask':15,'S_is_global_maximum':15 in maximum,'S_is_unique_maximum':maximum==[15],'eligible_subset_count':len(eligible),'eligible_hyperedge_connected':len(components)==1,'component_count':len(components),'connectivity_uses_all_eligible_hyperedges':True,'total_mass':str(totalmass),'whole_source_span_squared':str(spanall*spanall),'nonterminal_gap_span2_minus_35over33total':str(nonterminalgap),'nonterminal':nonterminalgap>0},'full_inner_old_involving_rectangles':fi,'full_outer_old_involving_rectangles':fo,'complete_residual_inner_rectangles':ri,'complete_residual_outer_rectangles':ro,'pureoutside_full_residual_cancellation_exact':True,'C_bounds':list(map(lambda z:str(F(z,L)),C)),'all_old_receiver_outer_rectangles_contained_in_C':True,'actual_full_plane_deletion_area':whole,'corner_regions':corners,'actual_noncorner_deletion_area':noncorner,'noncorner_definition':'C minus the union of four disjoint-interior corner squares; actual deletion outside C is empty. All boundaries have zero area.','integer_lower_and_upper_endpoints_add_exactly_in_decomposition':True,'SW_internal_local_mass':'1/4','SW_auxiliary_half_mass':str(halfP),'SW_strictly_exceeds_half_mass':F(sw['lower'])>halfP,'SW_excess_lower':str(F(sw['lower'])-halfP),'SW_main_corner_threshold':str(main),'SW_strictly_exceeds_main_corner_threshold':F(sw['lower'])>main,'whole_deletion_strict_below_3point5':F(whole['upper'])<F(7,2),'whole_3point5_minus_upper':str(F(7,2)-F(whole['upper'])),'algorithm_A':'vertical x-event direct y-interval union difference','algorithm_B':'compressed Cartesian midpoint-cell exact integer coverage','scope':'One fixed actual six-source measure retaining all four old positions. Local half-mass and main-corner tests and global deletion are separate. No grid collapse, parameter search or general C2 improvement. Qualification uses complete63subsets and actual bbox closures, including hyperedges for connectivity.','seconds':time.monotonic()-start}
 Path(args.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['pid','distinct_bbox_closure_count','qualification','actual_full_plane_deletion_area','corner_regions','actual_noncorner_deletion_area','SW_strictly_exceeds_half_mass','SW_strictly_exceeds_main_corner_threshold','whole_deletion_strict_below_3point5','seconds']},indent=2))
if __name__=='__main__':main()
