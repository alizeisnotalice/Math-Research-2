#!/usr/bin/env python3
"""R59 new weighted boundary kernel: independent half-axis bbox and source checks."""
import argparse,hashlib,importlib.util,json,math,os,sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True

def ceilq(x):return -(-x.numerator//x.denominator)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--kernel',type=Path,required=True);ap.add_argument('--axis-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ar=ap.parse_args();assert os.path.ismount('/Volumes/data')
 kh=sha(ar.kernel);sys.path.insert(0,str(ar.kernel.parent));sp=importlib.util.spec_from_file_location('r059_review_target',ar.kernel);target=importlib.util.module_from_spec(sp);sp.loader.exec_module(target)
 data=json.loads(ar.input.read_text());raw=[tuple(map(F,p)) for p in data['pts']];px,py=F(data['Px']),F(data['Py']);X=max(x for x,y,w in raw);axisrows=[];sourcerows=[]
 for side in ('left','right'):
  def setup(t):
   pts=[(X-x if side=='left' else x,y,w,(t if side=='left' else 1-t) if k in (0,51) else (F(1) if side=='left' else F(0))) for k,(x,y,w) in enumerate(raw)]
   anchor=X-3*px if side=='left' else -3*px;xmax=max(x if z>0 else x-px for x,y,w,z in pts);return pts,anchor,xmax
  def catalog(pts,anchor,C):
   atoms={};ids=set();R=F(15)
   for mirrored in (0,1):
    lower=C-R if mirrored else anchor-R;upper=2*C-anchor+R if mirrored else C+R
    for label,(a,b,w,z) in enumerate(pts):
     # Independent inclusive integer enumeration, then exact physical window filtering.
     for i in range(ceilq((lower-a)/px),min(0,(upper-a)//px)+1):
      ww=w*(z if i==0 else 1)
      if ww==0:continue
      xx=a+i*px;xx=2*C-xx if mirrored else xx
      for j in range(ceilq((-R-b)/py),(py+R-b)//py+1):
       yy=b+j*py;atoms[xx,yy]=atoms.get((xx,yy),F(0))+ww;ids.add((bool(mirrored),label,i,j))
   return atoms,ids
  for t in (F(0),F(1,2),F(1)):
   pts,anchor,xmax=setup(t);C=xmax+1
   got=target.setup(data,side,t);assert (got[3],got[4],got[5])==(pts,anchor,xmax)
   a,ids=catalog(pts,anchor,C);ga,gi=target.source_atoms(pts,px,py,anchor,C,F(15));assert ga==a and set(gi)==ids
   sourcerows.append(dict(side=side,t=str(t),probe_axis=str(C),atoms=len(a),source_dictionary_equals_target=True))
  pts,anchor,xmax=setup(F(1,2));p=ar.axis_dir/('r059-half-'+side+'-axis.json');ax=json.loads(p.read_text());assert ax['script_sha256']==kh and ax['input_sha256']==sha(ar.input)
  assert F(ax['actual_source_xmax'])==xmax;R=F(10)
  for bs in ax['bootstrap']:
   assert F(bs['R'])==R;nr=ceilq(2*R/py)+1;row=[]
   for label,(x,y,w,z) in enumerate(pts):
    for i in range((xmax-R-x)//px+1,1):
     ww=w*(z if i==0 else 1)
     if ww>0:row.append((i,label,x+i*px,y,ww))
   B=nr*sum(w for i,k,x,y,w in row);new=F(bs['next_strict_radius']);assert B==F(bs['mass_upper']) and nr==bs['rows'] and 4*new*new>B
   assert [[i,k,str(w)] for i,k,x,y,w in row]==bs['row_labels']
   if new>=R:break
   R=new
  assert R==F(ax['final_owner_radius']);atoms={}
  for i,k,x,y,w in row:
   for j in range(nr):atoms[x,y+j*py]=atoms.get((x,y+j*py),F(0))+w
  xs=sorted({x for x,y in atoms});ys=sorted({y for x,y in atoms});nx,ny=len(xs),len(ys);W=math.lcm(*(w.denominator for w in atoms.values()));xi={x:i+1 for i,x in enumerate(xs)};yi={y:j+1 for j,y in enumerate(ys)};pre=[[0]*(ny+1) for _ in range(nx+1)]
  for (x,y),w in atoms.items():pre[xi[x]][yi[y]]+=int(w*W)
  for i in range(1,nx+1):
   s=0
   for j in range(1,ny+1):s+=pre[i][j];pre[i][j]=pre[i-1][j]+s
  m0=max(w for (x,y),w in atoms.items() if x==xmax);q=m0/4;ct=dict(boxes=0,zero=0,empty=0,active=0,strict=0,equality=0,squared=0)
  for i,a in enumerate(xs):
   for k,b in enumerate(xs[i:],i):
    for j,c in enumerate(ys):
     for l,d in enumerate(ys[j:],j):
      ct['boxes']+=1;M=pre[k+1][l+1]-pre[i][l+1]-pre[k+1][j]+pre[i][j]
      if not M:ct['zero']+=1;continue
      if max((b-a)**2,(d-c)**2)>=F(M,W):ct['empty']+=1;continue
      ct['active']+=1;p0=F(M,4*W);g=xmax-a;L=p0-q-g*g
      if p0<=q:equal=p0==q and g==0
      elif L<=0:equal=False
      else:ct['squared']+=1;v=4*g*g*q-L*L;assert v>=0;equal=v==0
      ct['equality' if equal else 'strict']+=1
  lo,hi=F(ax['c_lower']),F(ax['c_upper']);assert 0<lo-xmax<hi-xmax and (lo-xmax)**2<=q<(hi-xmax)**2
  low,idslo=catalog(pts,anchor,lo);high,idshi=catalog(pts,anchor,hi);assert idshi<=idslo
  assert target.source_atoms(pts,px,py,anchor,lo,F(15))[0]==low and target.source_atoms(pts,px,py,anchor,hi,F(15))[0]==high
  cx=len({x for x,y in low});cy=len({y for x,y in low});K=cx*(cx+1)//2*cy*(cy+1)//2;err=2*(hi-lo)*py*K
  axisrows.append(dict(side=side,axis_sha256=sha(ar.axis_dir/('r059-half-'+side+'-axis.json')),catalog_atoms=len(atoms),nx=nx,ny=ny,counts=ct,exact_axis='actual_xmax+sqrt(q)',actual_xmax=str(xmax),q=str(q),author_interval_contains_exact_root=True,common_catalog=dict(atoms=len(low),nx=cx,ny=cy,all_bbox_count=K,high_labels_subset_low=True),true_axis_area_error=str(err)))
 assert sha(ar.kernel)==kh
 out=dict(status='PASS',python=sys.version,input_sha256=sha(ar.input),target_kernel_sha256=kh,script_sha256=sha(Path(__file__)),source_checks=sourcerows,independent_axis_checks=axisrows,scope='Half t only complete independent axes; setup/source dictionaries also checked at t0,1/2,1. No large halo area run; moving-anchor error proof separate.')
 ar.output.parent.mkdir(parents=True,exist_ok=True);ar.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
