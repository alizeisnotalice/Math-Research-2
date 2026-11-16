#!/usr/bin/env python3
"""R74 four controls x four strict source windows; exact regional checks.

For a corner receiver x and strict radius r<1/2, every captured source
lies in the appropriate OPEN P window. Maxblock caps exclude all strict
full or residual witnesses with r>=1/2. Deleting outside mass outside P
preserves caps, so full and residual halos are unchanged in that corner.
The geometric lemma gives equality; overlapping directed area brackets
below are diagnostic checks, not a substitute for its proof.
"""
import argparse,hashlib,importlib.util,json,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def interval(lo,hi,D):return {'lower':str(F(lo,D*D)),'upper':str(F(hi,D*D)),'lower_decimal':float(F(lo,D*D)),'upper_decimal':float(F(hi,D*D))}
def clipped(g,m,box):
 D=g['coordinate_scale'];bb=[]
 for x in box:assert (x*D).denominator==1;bb.append(int(x*D))
 job=next(j for j in g['all_maximum_block_evaluations'] if j['chosen_S']);fi=[r['inner'] for r in g['full_directed_rectangles']];fo=[r['outer'] for r in g['full_directed_rectangles']];ri=[r['inner'] for r in job['residual_directed_rectangles']];ro=[r['outer'] for r in job['residual_directed_rectangles']];lf,hr,ld=m.sweep(m.clip(fi,bb),m.clip(ro,bb));hf,lr,hd=m.sweep(m.clip(fo,bb),m.clip(ri,bb))
 return {'full_halo':interval(lf,hf,D),'residual_halo':interval(lr,hr,D),'deletion':interval(ld,hd,D)}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--archive',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();root=Path(a.archive);t=time.monotonic();enginepath=root/'computation/round072/actual_halo_geometry.py';regionpath=root/'computation/round073/regional_credit_diagnostic_v2.py';engine=module(enginepath,'frozen_r072_geometry');region=module(regionpath,'frozen_r073_clip');diagpath=root/'computation/results/round073-regional/four-control-regions-v2.json';old=json.loads(diagpath.read_text());w=F(24,25);e=F(1,25);L=(-F(1),F(1,2)-e);R=(F(1,2),w+1)
 windows={'SW':(*L,*L),'NW':(*L,*R),'SE':(*R,*L),'NE':(*R,*R)};corners={'SW':(F(-1,2),-e,F(-1,2),-e),'NW':(F(-1,2),-e,F(1),w+F(1,2)),'SE':(F(1),w+F(1,2),F(-1,2),-e),'NE':(F(1),w+F(1,2),F(1),w+F(1,2))};regionids={'SW':'x0y0','NW':'x0y2','SE':'x2y0','NE':'x2y2'};inputnames={'s25-only':'s25-only.json','r071-positive33':'r071-positive33.json','aligned225':'aligned225-actual-input.json','contracted-unique':'contracted-unique-actual-input.json'}
 cache={};cases=[];new_evaluations=0
 for c in old['cases']:
  name=c['name'];ip=root/'computation/inputs/r072'/inputnames[name];data=json.loads(ip.read_text());points=[tuple(map(F,p)) for p in data['positions']];masses=list(map(F,data['masses']));S=set(data['S_indices']);control={'name':name,'source_input_sha256':sha(ip),'corners':{}}
  for corner,P in windows.items():
   a0,a1,b0,b1=P;inside=lambda p:a0<p[0]<a1 and b0<p[1]<b1;onboundary=lambda p:a0<=p[0]<=a1 and b0<=p[1]<=b1 and (p[0] in (a0,a1) or p[1] in (b0,b1));keep=sorted(S)+[i for i in range(len(points)) if i not in S and inside(points[i])];removed=[i for i in range(len(points)) if i not in keep];boundary=[i for i in range(len(points)) if onboundary(points[i])];localS=[i for i in S if inside(points[i])]
   restricted={'threshold':'1','background':data['background'],'geometry':data['geometry'],'S_indices':list(range(len(S))),'positions':[data['positions'][i] for i in keep],'masses':[data['masses'][i] for i in keep]};raw=(json.dumps(restricted,indent=2)+'\n').encode();h=hashlib.sha256(raw).hexdigest();p=Path('computation/inputs/r074/localization')/(h+'.json')
   if h not in cache:
    p.write_bytes(raw)
    if len(removed)==0:
     gp=root/'computation/results/round072-actual'/f'{name}.json';g=json.loads(gp.read_text());origin='reused_identical_original_input'
    else:g=engine.evaluate(restricted,bits=72);origin='new_exact_restricted_input_geometry';new_evaluations+=1
    cache[h]=(g,origin)
   g,origin=cache[h];new=clipped(g,region,corners[corner]);previous=c['regions'][regionids[corner]];comparisons={}
   for metric in new:
    before=previous[metric];after=new[metric];lo=max(F(before['lower']),F(after['lower']));hi=min(F(before['upper']),F(after['upper']));assert lo<=hi;comparisons[metric]={'intervals_overlap':True,'overlap_lower':str(lo),'overlap_upper':str(hi),'maximum_endpoint_discrepancy':str(max(abs(F(before['lower'])-F(after['lower'])),abs(F(before['upper'])-F(after['upper']))))}
   control['corners'][corner]={'P_open':[str(x) for x in P],'receiver_corner':[str(x) for x in corners[corner]],'internal_sources_in_P':localS,'internal_mass_in_P':str(sum(masses[i] for i in localS)),'outside_original_count':len(points)-len(S),'outside_retained_count':len(keep)-len(S),'outside_removed_count':len(removed),'positive_sources_on_P_boundary_indices':boundary,'boundary_sources_excluded_from_outside_restriction':[i for i in boundary if i not in S],'restricted_input_path':str(p),'restricted_input_sha256':h,'new_positive_source_count':g['positive_source_count'],'restricted_input_keeps_S_maximum':g['chosen_S_is_maximum'],'geometry_origin':origin,'clipped_intervals':new,'comparisons_with_original_certified_corner':comparisons}
  cases.append(control)
 out={'schema':'r074-strict-corner-source-window-localization-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'frozen_geometry_source_path':str(enginepath),'frozen_geometry_source_sha256':sha(enginepath),'frozen_regional_source_path':str(regionpath),'frozen_regional_source_sha256':sha(regionpath),'source_original_region_path':str(diagpath),'source_original_region_sha256':sha(diagpath),'cases':cases,'distinct_restricted_inputs':len(cache),'new_complete_geometry_evaluations':new_evaluations,'seconds':time.monotonic()-t,'proof_interface':'All strict full/residual witnesses have r<1/2 by the shared square cap. At a SW corner point x in[-1/2,-e]^2, any captured source has each coordinate strictly between-1 and1/2-e. Reflection gives other windows, separated by a positive gap e. Removing outside mass outside its open window preserves all caps and neither adds nor removes a corner full/residual witness. AllS retained; only nu is clipped.','scope':'Four fixed controls, sixteen corner checks, exact directed intervals. Boundary source exclusion is explicit; all three regional measures are checked, not just their difference. New complete geometry is on restricted actual inputs only, with cache reuse; no search/optimization yet. Pointwise equality follows from analytic source-window lemma, not merely interval overlap.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':a.out,'pid':os.getpid(),'distinct_restricted_inputs':len(cache),'new_evaluations':new_evaluations,'seconds':out['seconds'],'all16corners_full_residual_deletion_consistent':True,'SWbaselines':[{'name':c['name'],'outside_retained':c['corners']['SW']['outside_retained_count'],'area':c['corners']['SW']['clipped_intervals']['deletion']} for c in cases]},indent=2))

if __name__=='__main__':main()
