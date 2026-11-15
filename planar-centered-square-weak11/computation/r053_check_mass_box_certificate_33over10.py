#!/usr/bin/env python3
"""Standalone exact replay; no heuristic/kernel/numpy/scipy dependency."""
from fractions import Fraction as Q
from pathlib import Path
import gzip,hashlib,json,math,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1];DIR=ROOT/'computation/results/r053-mass-adversary'
PATH=DIR/'all-mass-box-certificate-33over10.json.gz'

def main():
 assert os.path.ismount('/Volumes/data');start=time.perf_counter();d=json.loads(gzip.decompress(PATH.read_bytes()));nodes=d['nodes'];scale=d['coordinate_denominator'];assert scale==2**d['bits'] and scale%4==0
 seen=set();leafs=[];maximum=0
 def visit(k,inputlo,inputhi):
  nonlocal maximum
  assert k not in seen;seen.add(k);n=nodes[k]
  oldlo=list(inputlo);oldhi=list(inputhi);sl=sum(oldlo);sh=sum(oldhi)
  lo=[max(oldlo[i],1-sum(oldhi[j] for j in range(9) if j!=i)) for i in range(9)]
  hi=[min(oldhi[i],1-sum(oldlo[j] for j in range(9) if j!=i)) for i in range(9)]
  empty=sum(lo)>1 or sum(hi)<1 or any(a>b for a,b in zip(lo,hi))
  if n['kind']=='empty':assert empty;return
  assert not empty;assert list(map(Q,n['lo']))==lo and list(map(Q,n['hi']))==hi
  if n['kind']=='split':
   j=n['coordinate'];cut=Q(n['cut']);assert 0<=j<9 and lo[j]<cut<hi[j]
   ha=hi.copy();ha[j]=cut;lb=lo.copy();lb[j]=cut
   visit(n['left'],lo,ha);visit(n['right'],lb,hi);return
  assert n['kind']=='leaf';p=n['order'];assert sorted(p)==list(range(9));rect=[];prefs=[]
  for pos,j in enumerate(p):
   bound=min(sum(hi[t] for t in p[:pos+1]),1-sum(lo[t] for t in p[pos+1:]));assert 0<=bound<=1;prefs.append(str(bound))
   x,y=divmod(j,3);cuts=[0,scale//4,3*scale//4,scale];r=n['outer_rectangles_integer'][pos]
   h=cuts[x]-r[0];assert h>=0 and r==[cuts[x]-h,cuts[x+1]+h,cuts[y]-h,cuts[y+1]+h]
   # Verify claimed outward radius directly by exact squared inequalities.
   assert 4*h*h*bound.denominator>=scale*scale*bound.numerator
   if h:assert 4*(h-1)**2*bound.denominator<scale*scale*bound.numerator
   rect.append(tuple(r))
  assert prefs==n['prefix_upper']
  # Independent rectangle-cell enumeration, not generator sweep accumulation.
  xs=sorted(set(v for r in rect for v in r[:2]));ys=sorted(set(v for r in rect for v in r[2:]));num=0
  for a,b in zip(xs,xs[1:]):
   for c,e in zip(ys,ys[1:]):
    if any(x0<=a and b<=x1 and y0<=c and e<=y1 for x0,x1,y0,y1 in rect):num+=(b-a)*(e-c)
  assert num==n['area_upper_numerator'];assert n['area_denominator']==scale*scale;assert 10*num<33*scale*scale
  maximum=max(maximum,num);leafs.append(k)
 visit(0,[Q(0)]*9,[Q(1)]*9)
 assert len(seen)==len(nodes)==d['node_count'];assert len(leafs)==d['leaf_count'];assert str(Q(maximum,scale*scale))==d['uniform_area_upper_exact']
 result={'status':'PASS exact full mass-simplex BSP coverage and every leaf outward union bound','certificate_sha256':hashlib.sha256(PATH.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'nodes_verified':len(seen),'leaves_verified':len(leafs),'uniform_area_upper_exact':str(Q(maximum,scale*scale)),'strict_margin':str(Q(7,2)-Q(maximum,scale*scale)),'strict_margin_below_target':str(Q(33,10)-Q(maximum,scale*scale)),'method':'Exact fractions, independently checked split coverage and radius-square inequalities; independent full endpoint-grid area enumeration.','python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.perf_counter()-start}
 (DIR/'author-standalone-check-33over10.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
