#!/usr/bin/env python3
"""Exact finite BSP certificate: universal nine-cell ordered proxy <=7/2."""
from fractions import Fraction as F
from pathlib import Path
import ctypes,gzip,hashlib,io,json,math,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'computation/results/r053-mass-adversary'
LIB='/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r053_cell_union.dylib'
BITS=40;D=1<<BITS;TARGET=F(33,10);NODE_CAP=100000;BANDS=((0,F(1,4)),(F(1,4),F(3,4)),(F(3,4),1))
CELLS=[(BANDS[i//3][0],BANDS[i//3][1],BANDS[i%3][0],BANDS[i%3][1]) for i in range(9)]

def tighten(lo,hi):
 sl,sh=sum(lo),sum(hi)
 return [max(l,1-sh+h) for l,h in zip(lo,hi)],[min(h,1-sl+l) for l,h in zip(lo,hi)]
def union(rect):
 xs=sorted({x for r in rect for x in r[:2]});total=0
 for a,b in zip(xs,xs[1:]):
  ys=sorted((c,d) for x,y,c,d in rect if x<=a and b<=y);length=0;end=None
  for c,d in ys:
   if end is None or c>end:length+=d-c;end=d
   elif d>end:length+=d-end;end=d
  total+=(b-a)*length
 return total

def outer(lo,hi,order):
 s=F(0);remaining=sum(lo);rs=[];prefix=[]
 for j in order:
  s+=hi[j];remaining-=lo[j];u=min(s,1-remaining);assert 0<=u<=1
  # ceil(D sqrt(u)/2), certified by integer squares.
  q=u*D*D/4;n=math.isqrt(q.numerator//q.denominator)
  if n*n*q.denominator<q.numerator:n+=1
  assert n*n*q.denominator>=q.numerator and (n==0 or (n-1)**2*q.denominator<q.numerator)
  prefix.append(str(u));x0,x1,y0,y1=CELLS[j]
  # Retain every cell, including potentially empty ones: conservative.
  rs.append((int(x0*D)-n,int(x1*D)+n,int(y0*D)-n,int(y1*D)+n))
 return union(rs),rs,prefix

def main():
 assert os.path.ismount('/Volumes/data');start=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
 L=ctypes.CDLL(LIB);L.heuristic.restype=ctypes.c_double;DA=ctypes.c_double*9;IA=ctypes.c_int*9
 tree=[];leaves=[]
 def rec(lo,hi,depth):
  idx=len(tree);assert idx<NODE_CAP;tree.append(None);lo,hi=tighten(lo,hi)
  if sum(lo)>1 or sum(hi)<1 or any(l>h for l,h in zip(lo,hi)):
   tree[idx]={'kind':'empty'};return idx
  mid=[(l+h)/2 for l,h in zip(lo,hi)];sm=sum(mid);p=IA();L.heuristic(DA(*(float(x/sm) for x in mid)),p,3,53001);order=list(p)
  n,rs,prefix=outer(lo,hi,order)
  common={'lo':list(map(str,lo)),'hi':list(map(str,hi))}
  if TARGET.denominator*n<TARGET.numerator*D*D:
   leaf={'kind':'leaf',**common,'order':order,'prefix_upper':prefix,'area_upper_numerator':n,'area_denominator':D*D,'outer_rectangles_integer':rs,'depth':depth}
   tree[idx]=leaf;leaves.append(idx);return idx
  j=max(range(9),key=lambda j:hi[j]-lo[j]);m=(lo[j]+hi[j])/2;assert lo[j]<m<hi[j]
  h1=hi.copy();h1[j]=m;l2=lo.copy();l2[j]=m
  left=rec(lo,h1,depth+1);right=rec(l2,hi,depth+1)
  tree[idx]={'kind':'split',**common,'coordinate':j,'cut':str(m),'left':left,'right':right};return idx
 rec([F(0)]*9,[F(1)]*9,0)
 maxleaf=max(leaves,key=lambda k:tree[k]['area_upper_numerator']);maxnum=tree[maxleaf]['area_upper_numerator']
 result={'theorem':'For all w>=0 with sum(w)=1, some order pi of nine closed quarter/middle/quarter cells has union area <33/10, even retaining all zero-mass cells.','target':str(TARGET),'node_cap':NODE_CAP,'method':'Exact rational box tightening and binary splits; integer isqrt outward radius bounds; exact integer union area. Numerical heuristic only chooses a candidate order; each accepted bound is exact.','root_domain':'[0,1]^9 intersect sum(w)=1','bits':BITS,'coordinate_denominator':D,'nodes':tree,'leaf_count':len(leaves),'node_count':len(tree),'max_depth':max(tree[k]['depth'] for k in leaves),'maximum_leaf':maxleaf,'uniform_area_upper_exact':str(F(maxnum,D*D)),'uniform_area_upper_float':float(F(maxnum,D*D)),'strict_margin_below_7_2':str(F(7,2)-F(maxnum,D*D)),'strict_margin_below_target':str(TARGET-F(maxnum,D*D)),'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.perf_counter()-start,'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'order_heuristic_source_sha256':hashlib.sha256(Path(__file__).with_name('r053_cell_union.cpp').read_bytes()).hexdigest()}
 path=OUT/'all-mass-box-certificate-33over10.json.gz'
 with path.open('wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
   gz.write((json.dumps(result,indent=2)+'\n').encode())
 print(json.dumps({k:v for k,v in result.items() if k!='nodes'},indent=2));print('certificate_sha256',hashlib.sha256(path.read_bytes()).hexdigest())
if __name__=='__main__':main()
