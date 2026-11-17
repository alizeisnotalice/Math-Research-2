#!/usr/bin/env python3
"""Independent 9! enumeration for the r030 (L,d,e2) table."""
from __future__ import annotations
import itertools,json,sys
from collections import defaultdict
from pathlib import Path
OUT=Path(sys.argv[1]); REF=Path(sys.argv[2])
def main():
 d=json.loads(OUT.read_text());assert d['n']==3 and d['u128_checked']
 adj=[]
 for q in range(9):
  x,y=q%3,q//3;adj.append([Y*3+X for X,Y in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)) if 0<=X<3 and 0<=Y<3])
 got=defaultdict(int)
 for word in itertools.permutations(range(9)):
  pos=[0]*9
  for i,x in enumerate(word):pos[x]=i
  peaks=[x for x in range(9) if all(pos[x]<pos[y] for y in adj[x])]
  if len(peaks)!=2:continue
  v=word[0];p=peaks[0] if peaks[1]==v else peaks[1]
  z=[x for x in word if x!=v];z.remove(p);z=[p]+z
  L=next(i for i,x in enumerate(z,1) if x in adj[p])-1
  dx=abs(v%3-p%3);dy=abs(v//3-p//3);got[L,dx+dy,dx*dx+dy*dy]+=1
 tab=defaultdict(int)
 for h in d['holes']:
  for L,dist,e2,c in h['joint_H_L_d_e2']:tab[L,dist,e2]+=L*int(c)
 assert dict(got)==dict(tab),(len(got),len(tab))
 ref=json.loads(REF.read_text());dh=defaultdict(int);e=0
 for (L,dist,e2),c in got.items():dh[dist]+=c;e+=e2*c
 assert sum(got.values())==int(ref['A2']) and e==int(ref['squared_euclidean_total'])
 assert {str(k):v for k,v in dh.items()}=={k:int(v) for k,v in ref['distance_counts'].items()}
 out={'status':'passed','A2':str(sum(got.values())),'distance_counts':dict(dh),'e2_total':str(e),'joint_cells':len(got)}
 (OUT.parent/'n3-joint-audit.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
