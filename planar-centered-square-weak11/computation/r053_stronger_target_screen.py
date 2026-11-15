#!/usr/bin/env python3
import ctypes,time,json,os,hashlib
from pathlib import Path
assert os.path.ismount('/Volumes/data')
L=ctypes.CDLL('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/r053_cell_union.dylib');D=ctypes.c_double*9;I=ctypes.c_int*9
for n in ['upper','heuristic']:getattr(L,n).restype=ctypes.c_double
records=[]
for target in [3.4,3.35,3.3]:
 p=I();q=[(0,[0.]*9,[1.]*9)];total=closed=empty=0;start=time.perf_counter();maxdepth=0;maxu=0
 while q and total<100000:
  depth,lo,hi=q.pop();total+=1;sl=sum(lo);sh=sum(hi)
  if sl>1+1e-14 or sh<1-1e-14:empty+=1;continue
  ol=lo;oh=hi;lo=[max(l,1-sh+h) for l,h in zip(ol,oh)];hi=[min(h,1-sl+l) for l,h in zip(ol,oh)]
  if any(l>h+1e-14 for l,h in zip(lo,hi)):empty+=1;continue
  mid=[(l+h)/2 for l,h in zip(lo,hi)];sm=sum(mid);w=[m/sm for m in mid]
  L.heuristic(D(*w),p,3,53001);u=L.upper(D(*lo),D(*hi),p)
  if u<target-1e-12:closed+=1;maxu=max(maxu,u);continue
  j=max(range(9),key=lambda j:hi[j]-lo[j]);m=(hi[j]+lo[j])/2;h1=hi.copy();h1[j]=m;l2=lo.copy();l2[j]=m;q.append((depth+1,lo,h1));q.append((depth+1,l2,hi));maxdepth=max(depth,maxdepth)
 records.append({'scope':'floating feasibility screen only, full simplex, not exact certificate','target':target,'node_cap':100000,'total':total,'closed':closed,'empty':empty,'open':len(q),'maxdepth':maxdepth,'elapsed':time.perf_counter()-start,'max_closed_bound':maxu});print(records[-1],flush=True)
 if q:break
out=Path(__file__).resolve().parent/'results/r053-mass-adversary/stronger-target-screen.json';out.write_text(json.dumps({'records':records,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
