#!/usr/bin/env python3
"""Three specified inputs; mixed full labels only, complete residual halo.

Radius intervals are exact rational directed 72-bit enclosures. Two frozen
integer union algorithms independently evaluate each directed area endpoint.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib.util, json, math, os, platform, sys, time
sys.dont_write_bytecode = True
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
ROOT=Path(__file__).resolve().parents[2]
CORE=ROOT/'computation/round075/certify_three_atom_collapse.py'
assert sha(CORE)=='a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456'
spec=importlib.util.spec_from_file_location('frozen_r075_area',CORE)
area=importlib.util.module_from_spec(spec);spec.loader.exec_module(area)
def run_case(data,bits):
 p=[tuple(map(F,x)) for x in data['positions']];m=list(map(F,data['masses']))
 finite=data['clip']!=['SW'];clip=list(map(F,data['clip'])) if finite else None
 den=[v.denominator for q in p for v in q]+[v.denominator for v in m]+([] if clip is None else [v.denominator for v in clip])
 scale=math.lcm(*den)*(1<<bits);rows=[]
 for mask in range(1,8):
  jj=[j for j in range(3) if mask>>j&1];mass=sum(m[j] for j in jj)
  box=(min(p[j][0] for j in jj),max(p[j][0] for j in jj),min(p[j][1] for j in jj),max(p[j][1] for j in jj))
  span=max(box[1]-box[0],box[3]-box[2]);row={'mask':mask,'mass':str(mass),'bbox':list(map(str,box)),'eligible':mass>span*span,'span_squared':str(span*span)}
  if row['eligible']:
   n=mass*scale*scale;assert n.denominator==1
   k=math.isqrt(n.numerator//4);up=k if 4*k*k==n else k+1
   assert 4*k*k<=n<=4*up*up
   aa,bb,cc,dd=[int(v*scale) for v in box]
   row.update(radius_lower_units=k,radius_upper_units=up,inner=[bb-k,aa+k,dd-k,cc+k],outer=[bb-up,aa+up,dd-up,cc+up])
  rows.append(row)
 def clipped(rr):
  if not finite:return area.SW(rr)
  a,b,c,d=[int(v*scale) for v in clip];out=[]
  for u,v,s,t in rr:
   z=[max(a,u),min(b,v),max(c,s),min(d,t)]
   if z[0]<z[1] and z[2]<z[3]:out.append(z)
  return out
 def labels(key,full):return clipped([r[key] for r in rows if r['eligible'] and r['mask'] in ((3,5,7) if full else (2,4,6))])
 fi,fo,ri,ro=labels('inner',True),labels('outer',True),labels('inner',False),labels('outer',False)
 lo=area.event_area(fi,ro);hi=area.event_area(fo,ri)
 lo2,nlo=area.midpoint_area(fi,ro);hi2,nhi=area.midpoint_area(fo,ri)
 assert (lo,hi)==(lo2,hi2) and lo<=hi
 low,high=F(lo,scale*scale),F(hi,scale*scale);target=F(data['target'])
 return {'id':data['id'],'scope':data['scope'],'scale':scale,'bits':bits,'full_masks':[3,5,7],'residual_masks':[2,4,6],'old_singleton_included':False,'catalog':rows,'full_inner':fi,'full_outer':fo,'residual_inner':ri,'residual_outer':ro,'area':{'lower':str(low),'upper':str(high),'lower_decimal':float(low),'upper_decimal':float(high),'width':str(high-low)},'target':str(target),'strict_above':low>target,'strict_below':high<target,'excess_lower':str(low-target),'two_algorithms_identical':True,'midpoint_cells':[nlo,nhi]}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();t=time.monotonic();ip=Path(a.input);data=json.loads(ip.read_text());out={'schema':'r077-mixed-proxy-directed-area-v1','pid':os.getpid(),'python':sys.version,'platform':platform.platform(),'script_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'frozen_area_core_sha256':sha(CORE),'algorithm_A':'vertical x-event/direct y-interval exact union difference','algorithm_B':'Cartesian compressed midpoint exact coverage','cases':[run_case(c,data['bits']) for c in data['cases']]}
 assert out['cases'][0]['strict_above'] and out['cases'][1]['strict_below'] and out['cases'][2]['strict_above']
 out['seconds']=time.monotonic()-t;Path(a.out).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({c['id']:c['area']|{'above':c['strict_above'],'below':c['strict_below']} for c in out['cases']},indent=2))
if __name__=='__main__':main()
