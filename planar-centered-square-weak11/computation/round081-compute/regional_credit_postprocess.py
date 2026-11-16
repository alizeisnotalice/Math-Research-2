#!/usr/bin/env python3
"""Reclip saved endpoint geometry only; no new inputs, caps or source search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,json,os,sys,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
CORE=ROOT/'computation/round075/certify_three_atom_collapse.py';assert sha(CORE)=='a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456'
spec=importlib.util.spec_from_file_location('frozen_r075',CORE);area=importlib.util.module_from_spec(spec);spec.loader.exec_module(area)
def main():
 started=time.monotonic();rp=ROOT/'computation/results/round081-compute';prod=rp/'production';summary=json.loads((prod/'summary.json').read_text());out=[]
 for endpoint in summary['endpoints']:
  ip=prod/endpoint['result'];assert sha(ip)==endpoint['result_sha256'];data=json.loads(ip.read_text());L=data['integer_scale'];full=data['full_bbox_catalog'];res=data['residual_bbox_catalog']
  def rr(cat,key):return [r[key] for r in cat if r['eligible']]
  def clip(rr,box):
   out=[]
   for a,b,c,d in rr:
    z=[max(a,box[0]),min(b,box[1]),max(c,box[2]),min(d,box[3])]
    if z[0]<z[1] and z[2]<z[3]:out.append(z)
   return out
  corners=[list(map(lambda x:int(F(x)*L),data['five_regions'][j]['bounds'])) for j in ['SW','SE','NE','NW']]
  def union_bracket(cat,box,excluded):
   FI=clip(rr(cat,'inner'),box);FO=clip(rr(cat,'outer'),box);lo=area.event_area(FI,excluded);hi=area.event_area(FO,excluded);lo2,_=area.midpoint_area(FI,excluded);hi2,_=area.midpoint_area(FO,excluded);assert (lo,hi)==(lo2,hi2)
   return [F(lo,L*L),F(hi,L*L)]
  regions={};quota=F(data['aggregate_corner_quota']);excesslo=excesshi=reservelo=reservehi=F(0)
  for name,reg in data['five_regions'].items():
   box=list(map(lambda x:int(F(x)*L),reg['bounds']));excluded=corners if name=='noncorner' else []
   size=F((box[1]-box[0])*(box[3]-box[2]),L*L)-(sum((F((q[1]-q[0])*(q[3]-q[2]),L*L) for q in corners),F(0)) if name=='noncorner' else F(0))
   HF=union_bracket(full,box,excluded);R=union_bracket(res,box,excluded);U=[size-HF[1],size-HF[0]];D=[F(reg['actual_deletion'][k]) for k in ['lower','upper']];deficit=[size-D[1],size-D[0]]
   assert max(deficit[0],U[0]+R[0])<=min(deficit[1],U[1]+R[1])
   value={'region_area':str(size),'full_halo_area_interval':list(map(str,HF)),'complete_residual_area_interval':list(map(str,R)),'uncovered_full_area_interval':list(map(str,U)),'actual_deletion_deficit_interval':list(map(str,deficit)),'actual_identity':'region minus D = uncovered-full disjoint-union residual','deficit_interval_overlaps_U_plus_R_interval':True}
   if name!='noncorner':
    ex=[max(F(0),x-quota/4) for x in D];re=[max(F(0),quota/4-D[1]),max(F(0),quota/4-D[0])];excesslo+=ex[0];excesshi+=ex[1];reservelo+=re[0];reservehi+=re[1];value.update(local_quota_excess_interval=list(map(str,ex)),local_quota_reserve_interval=list(map(str,re)))
   regions[name]=value
  out.append({'case_id':data['case_id'],'t':data['t'],'endpoint_result_sha256':sha(ip),'regions':regions,'sum_corner_excess_interval':list(map(str,[excesslo,excesshi])),'sum_corner_reserve_interval':list(map(str,[reservelo,reservehi])),'reserve_minus_excess_interval':list(map(str,[reservelo-excesshi,reservehi-excesslo])),'reserve_minus_excess_decimal_interval':[float(reservelo-excesshi),float(reservehi-excesslo)],'sum_corner_excess_decimal_interval':[float(excesslo),float(excesshi)],'sum_corner_reserve_decimal_interval':[float(reservelo),float(reservehi)]})
 result={'schema':'r081-five-region-shared-credit-v1','pid':os.getpid(),'python':sys.version,'source_sha256':sha(__file__),'production_summary_sha256':sha(prod/'summary.json'),'cases':out,'scope':'Only reclip the two saved endpoint geometries. Complete full/residual halos and true deletion retained. Local excess/reserve is bookkeeping, not an independent four-corner theorem. No parameter exploration or uniform B&B.','seconds':time.monotonic()-started};op=rp/'regional-credit-diagnostics.json';assert not op.exists();op.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'cases':[{k:c[k] for k in ['case_id','t','sum_corner_excess_decimal_interval','sum_corner_reserve_decimal_interval','reserve_minus_excess_decimal_interval']} for c in out],'seconds':result['seconds']},indent=2))
if __name__=='__main__':main()
