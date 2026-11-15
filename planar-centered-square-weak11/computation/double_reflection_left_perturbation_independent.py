#!/usr/bin/env python3
"""R54 independent seam-only low-axis catalogue and perturbation replay."""
import argparse,hashlib,json,os,sys
from pathlib import Path
from fractions import Fraction as F

def ceilq(x):return -(-x.numerator//x.denominator)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--axis',type=Path,required=True);ap.add_argument('--error',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert os.path.ismount('/Volumes/data')
 data=json.loads(a.input.read_text());ax=json.loads(a.axis.read_text());er=json.loads(a.error.read_text());px,py=F(data['Px']),F(data['Py']);pts=[tuple(map(F,p)) for p in data['pts']];lo,hi=F(ax['c_edge_lower']),F(ax['c_edge_upper']);delta=hi-lo
 assert lo>max(x for x,y,m in pts) and delta==F(1,2**64)
 assert er['input_sha256']==hashlib.sha256(a.input.read_bytes()).hexdigest() and er['seam_sha256']==hashlib.sha256(a.axis.read_bytes()).hexdigest()
 by={r['state']:r for r in er['states']};rows=[];errors=[]
 for state,R in [('seam-old',F(10)),('seam-new',F(15))]:
  def catalog(c):
   atoms={};labels=set()
   for mirrored in range(2 if state=='seam-new' else 1):
    lower=c-R if mirrored else c-15-R
    upper=c+15+R if mirrored else c+R
    for k,(x,y,m) in enumerate(pts):
     for i in range(ceilq((lower-x)/px),min(0,(upper-x)//px)+1):
      for j in range(ceilq((-R-y)/py),(py+R-y)//py+1):
       xx=x+i*px;xx=2*c-xx if mirrored else xx;yy=y+j*py
       atoms[xx,yy]=atoms.get((xx,yy),F(0))+m;labels.add((mirrored,k,i,j))
   return atoms,labels
  atoms,lowids=catalog(lo);_,highids=catalog(hi);assert highids<=lowids
  nx=len({x for x,y in atoms});ny=len({y for x,y in atoms});K=nx*(nx+1)//2*ny*(ny+1)//2;err=2*delta*py*K
  r=by[state];assert len(atoms)==r['atoms'] and nx==r['nx'] and ny==r['ny'] and K==r['all_bbox_upper_count'] and err==F(r['area_change_upper'])
  rows.append(dict(state=state,atoms=len(atoms),nx=nx,ny=ny,complete_bbox_count=K,area_change_upper=str(err),high_axis_labels_subset_low=True));errors.append(err)
 assert sum(errors)==F(er['gamma_change_upper'])
 result=dict(status='PASS',python=sys.version,arithmetic='Fraction; exact inclusive lattice ranges and full bbox count',rows=rows,gamma_change_upper=str(sum(errors)),gamma_change_upper_display=float(sum(errors)),source_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),axis_sha256=hashlib.sha256(a.axis.read_bytes()).hexdigest(),author_error_sha256=hashlib.sha256(a.error.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Only seam old/new perturbation required by vertically periodic R54; no corner claims and no large halo-area rerun')
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
