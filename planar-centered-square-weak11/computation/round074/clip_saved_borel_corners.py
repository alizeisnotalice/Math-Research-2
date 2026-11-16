#!/usr/bin/env python3
"""R74 single prediagnostic: clip frozen R72 Borel leaf certificates.

No classifications, distance events, square roots or subdivisions are run.
All area numerators are integers in the inherited coordinate units.
"""
import argparse,gzip,hashlib,json,os,platform,sys,time
from pathlib import Path
from fractions import Fraction as F

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(lo,hi,D):
 lo=F(lo,D*D);hi=F(hi,D*D);assert 0<=lo<=hi
 return {'lower':str(lo),'upper':str(hi),'lower_decimal':float(lo),'upper_decimal':float(hi),'width':str(hi-lo)}
def run(root,n):
 t=time.monotonic();rp=root/f'computation/results/round072-borel/n{n}.json';cp=root/f'computation/results/round072-borel/n{n}.cells.json.gz';ap=root/f'computation/results/round072-independent/borel-n{n}.json';summary=json.loads(rp.read_text());audit=json.loads(ap.read_text());assert sha(cp)==summary['certificate_sha256']==audit['certificate_sha256'];assert sha(rp)==audit['summary_sha256']
 cert=json.loads(gzip.decompress(cp.read_bytes()));D=cert['integer_coordinate_scale'];assert D==summary['integer_coordinate_scale'];C=cert['receiver_domain_units'];jl,jh=cert['J_axis_units'];e=F(1,n*n);assert F(jl,D)==-e and F(jh,D)==1
 boxes={'left_bottom':(C[0],jl,C[2],jl),'left_top':(C[0],jl,jh,C[3]),'right_bottom':(jh,C[1],C[2],jl),'right_top':(jh,C[1],jh,C[3])};sums={k:{-1:0,0:0,1:0} for k in boxes};hits={k:{-1:0,0:0,1:0} for k in boxes};global_area={-1:0,0:0,1:0};counts={-1:0,0:0,1:0}
 for leaf in cert['leaves']:
  state=leaf['state'];assert state in global_area;a,b,c,d=leaf['cell_units'];assert a<b and c<d;global_area[state]+=(b-a)*(d-c);counts[state]+=1
  for key,(u,v,s,tt) in boxes.items():
   width=min(b,v)-max(a,u);height=min(d,tt)-max(c,s)
   if width>0 and height>0:sums[key][state]+=width*height;hits[key][state]+=1
 assert {str(k):v for k,v in counts.items()}==summary['leaf_counts'];assert sum(global_area.values())==(C[1]-C[0])*(C[3]-C[2]);assert F(global_area[1],D*D)==F(summary['area_lower']);assert F(global_area[1]+global_area[0],D*D)==F(summary['area_upper'])
 corners={}
 for key,box in boxes.items():
  area=(box[1]-box[0])*(box[3]-box[2]);assert sum(sums[key].values())==area;assert F(area,D*D)==(F(1,2)-e)**2
  corners[key]={'bounds':[str(F(a,D)) for a in box],'area':str(F(area,D*D)),'deletion':enc(sums[key][1],sums[key][1]+sums[key][0],D),'excluded_area':str(F(sums[key][-1],D*D)),'unresolved_area':str(F(sums[key][0],D*D)),'intersecting_leaf_counts_by_state':{str(s):v for s,v in hits[key].items()}}
 total={s:sum(sums[k][s] for k in boxes) for s in global_area};interval=enc(total[1],total[1]+total[0],D);threshold=F(1,2)+3*e*e;passed=F(interval['upper'])<=threshold;failed=F(interval['lower'])>threshold
 return {'n':n,'density_t':summary['density_t'],'e':str(e),'C':summary['receiver_domain'],'J':summary['J'],'coordinate_scale':D,'source_summary_path':str(rp),'source_summary_sha256':sha(rp),'source_leaf_certificate_path':str(cp),'source_leaf_certificate_sha256':sha(cp),'source_independent_audit_path':str(ap),'source_independent_audit_sha256':sha(ap),'inherited_audit_status':audit['status'],'inherited_classified_leaf_replay_counts':audit['independently_replayed_classified_leaves'],'inherited_audit_scope':audit['scope'],'source_total_leaf_counts':counts,'corners':corners,'total_corner_area':str(4*(F(1,2)-e)**2),'total_corner_deletion':interval,'corner_only_threshold_half_plus_3e_squared':str(threshold),'threshold_decimal':float(threshold),'upper_margin_below_threshold':str(threshold-F(interval['upper'])),'lower_excess_above_threshold':str(F(interval['lower'])-threshold),'decision':'PASSES_THIS_BOREL_CONTROL' if passed else 'FAILS_ON_THIS_BOREL_CONTROL' if failed else 'UNRESOLVED_BY_SAVED_LEAVES','new_leaf_classifications':0,'new_subdivisions':0,'seconds':time.monotonic()-t}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--archive',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();t=time.monotonic();cases=[run(Path(a.archive),n) for n in (5,17)];out={'schema':'r074-saved-borel-corner-clipping-v1','python':sys.version,'platform':platform.platform(),'pid':os.getpid(),'source_sha256':sha(Path(__file__)),'cases':cases,'seconds':time.monotonic()-t,'method':'For each frozen leaf, integer rectangle intersection with each corner. State1 intersection charged to lower and upper; state0 to upper only; state-1 to neither. Recomputed only clipped bookkeeping and verified old global sums, not old classifications or partitions.','scope':'Actual Borel outside density with inherited proof that residual halo is empty; thus corner Delta equals full halo there. No finite-atomic conclusion, area transfer, new classification, or general corner bound. n5 inherited all classified-leaf independent replay; n17 inherited only32 classified samples and full partition/bookkeeping review.'};Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'output':a.out,'pid':os.getpid(),'seconds':out['seconds'],'cases':[{'n':c['n'],'corners':{k:v['deletion'] for k,v in c['corners'].items()},'total':c['total_corner_deletion'],'threshold':c['corner_only_threshold_half_plus_3e_squared'],'decision':c['decision']} for c in cases]},indent=2))

if __name__=='__main__':main()
