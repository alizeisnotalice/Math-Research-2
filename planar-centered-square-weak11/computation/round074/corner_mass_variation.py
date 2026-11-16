#!/usr/bin/env python3
"""Three true-corner mass directions from one selected finite seed.

Finite differences/HiGHS only propose directions. Common-denominator rational
rounding and all global squarecaps decide admissibility. Full72bit complete
halo-minus-residual corner intervals decide every accepted improvement.
"""
import argparse,hashlib,importlib.util,json,math,os,platform,sys,time,warnings
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import numpy as np
import scipy
from scipy import sparse
from scipy.optimize import linprog
sys.dont_write_bytecode=True
U=200;D=40000
S=[(48*i,48*j) for i in range(5) for j in range(5)]
O=[(48*i,48*j) for i in range(-4,2) for j in range(-4,2) if not(0<=i<=1 and 0<=j<=1)]+[(-60,47)]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mod(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def ids(mask):
 out=[]
 while mask:b=mask&-mask;out.append(b.bit_length()-1);mask-=b
 return out
def bbox(ps):return (min(p[0] for p in ps),max(p[0] for p in ps),min(p[1] for p in ps),max(p[1] for p in ps))
def closure_masks(points):
 axes=[]
 for axis in (0,1):
  groups=defaultdict(int)
  for j,p in enumerate(points):groups[p[axis]]|=1<<j
  xx=sorted(groups);aa=[]
  for i in range(len(xx)):
   mask=0
   for x in xx[i:]:mask|=groups[x];aa.append(mask)
  axes.append(aa)
 return sorted({a&b for a in axes[0] for b in axes[1]}-{0})
def csr(rows,n):
 rr=[];cc=[]
 for i,row in enumerate(rows):rr.extend([i]*len(row));cc.extend(row)
 return sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(rows),n))
def build():
 n=len(O);maskO=(1<<n)-1;P=O+S;full=[];caps={}
 for mask in closure_masks(P):
  jj=ids(mask);box=bbox([P[j] for j in jj]);outmask=mask&maskO;rr=ids(outmask);ns=(mask>>n).bit_count();span=max(box[1]-box[0],box[3]-box[2]);full.append((rr,F(ns,25),box))
  if outmask:
   rhs=max(D,span*span)-1600*ns
   if outmask not in caps or rhs<caps[outmask][0]:caps[outmask]=(rhs,rr)
 caprows=[v[1] for v in caps.values()];capnum=[v[0] for v in caps.values()];resrows=caprows;resboxes=[bbox([O[j] for j in rr]) for rr in resrows]
 return {'caprows':caprows,'capnum':capnum,'A':csr(caprows,n),'b':np.array(capnum)/D,'FC':csr([a[0] for a in full],n),'sigma':np.array([float(a[1]) for a in full]),'Fbox':np.array([a[2] for a in full])/U,'RC':csr(resrows,n),'Rbox':np.array(resboxes)/U,'full_closures':len(full),'residual_closures':len(resrows)}
def float_rectangles(masses,boxes):
 span=np.maximum(boxes[:,1]-boxes[:,0],boxes[:,3]-boxes[:,2]);ok=masses>span*span;B=boxes[ok];r=np.sqrt(masses[ok])/2
 if not len(B):return []
 aa=np.maximum(B[:,1]-r,-.5);bb=np.minimum(B[:,0]+r,-.04);cc=np.maximum(B[:,3]-r,-.5);dd=np.minimum(B[:,2]+r,-.04);good=(aa<bb)&(cc<dd)
 return list(zip(aa[good],bb[good],cc[good],dd[good]))
def float_objective(g,y,region):return float(region.sweep(float_rectangles(g['FC']@y+g['sigma'],g['Fbox']),float_rectangles(g['RC']@y,g['Rbox']))[2])
def rationalize(y,g):
 q=[F(round(max(0.,float(a))*100000000),100000000) for a in y];scale=F(1)
 for row,b in zip(g['caprows'],g['capnum']):
  z=sum(q[j] for j in row)
  if z>F(b,D):scale=min(scale,F(b,D)/z)
 q=[a*scale for a in q];assert all(sum(q[j] for j in row)<=F(b,D) for row,b in zip(g['caprows'],g['capnum']));return q,scale
def input_data(q):
 pp=[(p,a) for p,a in zip(O,q) if a>0]
 return {'threshold':'1','background':'Lebesgue R2','geometry':'centered axis-parallel common-radius squares','S_indices':list(range(25)),'positions':[[str(F(x,U)),str(F(y,U))] for x,y in S]+[[str(F(x,U)),str(F(y,U))] for (x,y),a in pp],'masses':['1/25']*25+[str(a) for p,a in pp]}
def exact(q,engine,region):
 data=input_data(q);r=engine.evaluate(data,bits=72);assert r['maximum_mass']=='1' and r['chosen_S_is_maximum'];s=r['coordinate_scale'];box=[]
 for a in (F(-1,2),F(-1,25),F(-1,2),F(-1,25)):assert (a*s).denominator==1;box.append(int(a*s))
 job=next(a for a in r['all_maximum_block_evaluations'] if a['chosen_S']);fi=[a['inner'] for a in r['full_directed_rectangles']];fo=[a['outer'] for a in r['full_directed_rectangles']];ri=[a['inner'] for a in job['residual_directed_rectangles']];ro=[a['outer'] for a in job['residual_directed_rectangles']];lf,hr,ld=region.sweep(region.clip(fi,box),region.clip(ro,box));hf,lr,hd=region.sweep(region.clip(fo,box),region.clip(ri,box));enc=lambda a,b:engine.enc(a,b,s*s)
 return {'corner_Delta':enc(ld,hd),'corner_full_halo':enc(lf,hf),'corner_residual_halo':enc(lr,hr),'actual_state':{k:r[k] for k in ['positive_source_count','maximum_count','full_closure_count','full_eligible_closure_count','eligible_component_count','actual_total_mass','actual_span','nonterminal','nonterminal_gap_span_squared_minus_35total_over33']},'chosen_S_global_Delta':job['actual_Delta'],'full_global_halo':r['full_halo']}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--archive',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();start=time.monotonic();root=Path(a.archive);ep=root/'computation/round072/actual_halo_geometry.py';rp=root/'computation/round073/regional_credit_diagnostic_v2.py';engine=mod(ep,'engine');region=mod(rp,'region');g=build();q=[F(0)]*len(O);q[O.index((-96,48))]=q[O.index((48,-96))]=F(8,25);current=exact(q,engine,region);seed=current;ledger=[];accepted=[];eps=1e-5;trust=F(1,16)
 for step in range(3):
  y=np.array(list(map(float,q)));base=float_objective(g,y,region);gradient=[]
  for j in range(len(O)):
   plus=y.copy();plus[j]+=eps
   if y[j]>eps:
    minus=y.copy();minus[j]-=eps;gradient.append((float_objective(g,plus,region)-float_objective(g,minus,region))/(2*eps))
   else:gradient.append((float_objective(g,plus,region)-base)/eps)
  bounds=[(max(0.,x-float(trust)),min(1.,x+float(trust))) for x in y]
  with warnings.catch_warnings(record=True) as ww:lp=linprog(-np.array(gradient),A_ub=g['A'],b_ub=g['b'],bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'threads':1})
  entry={'direction':step,'float_baseline':base,'exact_baseline':current['corner_Delta'],'gradient':gradient,'LPstatus':lp.status,'LPmessage':lp.message,'warnings':[str(w.message) for w in ww],'line_trials':[]}
  if not lp.success:entry['stop_reason']='NO_LP_PROPOSAL';ledger.append(entry);break
  best=None
  for alpha in [F(1),F(1,2),F(1,4),F(1,8)]:
   proposal=y+float(alpha)*(lp.x-y);qq,scale=rationalize(proposal,g);evaluation=exact(qq,engine,region);improved=F(evaluation['corner_Delta']['lower'])>F(current['corner_Delta']['upper']);trial={'alpha':str(alpha),'weights':[str(v) for v in qq],'rationalization_scale':str(scale),'exact_evaluation':evaluation,'strict_improvement':improved};entry['line_trials'].append(trial)
   if improved and (best is None or F(evaluation['corner_Delta']['lower'])>F(best[1]['corner_Delta']['lower'])):best=(qq,evaluation,alpha)
  if best is None:entry['stop_reason']='NO_CERTIFIED_TRUE_AREA_IMPROVEMENT';ledger.append(entry);break
  q,current,alpha=best;accepted.append({'direction':step,'alpha':str(alpha),'corner_Delta':current['corner_Delta']});entry['accepted_alpha']=str(alpha);ledger.append(entry);print(json.dumps({'direction':step,'accepted_alpha':str(alpha),'corner_Delta':current['corner_Delta'],'positive_outside':sum(v>0 for v in q)}),flush=True)
 out={'schema':'r074-one-seed-true-corner-mass-variation-v1','python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'frozen_geometry_source_sha256':sha(ep),'frozen_regional_source_sha256':sha(rp),'parameters':{'n':5,'outside_position_count':len(O),'outside_positions_units_1over200':O,'full_global_S_source_count':25,'all_global_maxcap_rows':len(g['caprows']),'full_fixed_source_closures':g['full_closures'],'residual_fixed_source_closures':g['residual_closures'],'seed':'cross-pair,8/25 each','maximum_direction_calls':3,'actual_direction_calls':len(ledger),'finite_difference_step':eps,'trust_radius':str(trust),'line_steps':['1','1/2','1/4','1/8'],'mass_rounding_denominator':100000000,'acceptance':'new exact72bit lower>current exact72bit upper','seed_random':None,'LPthreads_requested':1},'seed_exact':seed,'ledger':ledger,'accepted':accepted,'final_weights':[str(v) for v in q],'final_exact':current,'strict_auxiliary_violation':F(current['corner_Delta']['lower'])>F(2,25),'strict_main_violation':F(current['corner_Delta']['lower'])>F(631,5000),'seconds':time.monotonic()-start,'scope':'Single chosen finite-position mass model. Floats propose directions only. All global maxcaps and complete true full/residual corner area checked for every rational line trial. No positions/n expansion, R72 proxy, local/global optimality or general inequality claim.'};p=Path(a.out);p.write_text(json.dumps(out,indent=2)+'\n');ip=Path('computation/inputs/r074/mass-variation/final-corner-input.json');ip.write_text(json.dumps(input_data(q),indent=2)+'\n');print(json.dumps({'output':str(p),'final_input':str(ip),'input_sha256':sha(ip),'seconds':out['seconds'],'final_corner':current['corner_Delta'],'auxiliary_violation':out['strict_auxiliary_violation'],'main_violation':out['strict_main_violation']},indent=2))
if __name__=='__main__':main()
