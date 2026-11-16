from fractions import Fraction as F
from pathlib import Path
import json,sys,time,platform,os,hashlib
cfg=json.loads(Path(sys.argv[1]).read_text());n=cfg['n'];c=F(cfg['density']);N=cfg['receiver_mesh'];start=time.time();e=F(1,n*n);h=F(n+1,n*n);w=1-e;p=F(1,2)-e;coords=[i*h for i in range(n)];q=n//2
length=lambda a,b,u,v:max(F(0),min(b,v)-max(a,u))
cuts=sorted({F(-1),w}|{v for a in coords+[F(-1),p,-e] for v in [a,a-1] if -1<=v<=w})
vertices=[];ratios=[]
for u in cuts:
 for v in cuts:
  k=sum(u<=s<=u+1 for s in coords);l=sum(v<=s<=v+1 for s in coords)
  A=length(u,u+1,-1,p)*length(v,v+1,-1,p)-length(u,u+1,-e,p)*length(v,v+1,-e,p)
  slack=1-F(k*l,n*n)-c*A;assert slack>=0,(u,v,k,l,A,slack)
  if A:ratios.append(((1-F(k*l,n*n))/A,u,v))
  vertices.append({'origin':[str(u),str(v)],'internal_count':k*l,'outside_area':str(A),'slack':str(slack)})
cap,cu,cv=min(ratios);assert cap==c
side=F(1,2)-e;step=side/N;cells=[];tests=0
for i in range(N):
 x0=-F(1,2)+i*step;x1=x0+step
 for j in range(N):
  y0=-F(1,2)+j*step;y1=y0+step
  for r in sorted(set([s-x0 for s in coords[:q]]+[s-y0 for s in coords[:q]])):
   if r>=F(1,2):continue
   k=sum(s<=x0+r for s in coords);l=sum(s<=y0+r for s in coords)
   if not k*l:continue
   Jmax=(x1+r+e)*(y1+r+e)
   margin=F(k*l,n*n)-4*(1-c)*r*r-c*Jmax;tests+=1
   if margin>0:
    assert x0-r>-1 and y0-r>-1 and x1+r<p and y1+r<p
    cells.append({'i':i,'j':j,'closed_radius':str(r),'guaranteed_internal_count':k*l,'strict_margin':str(margin)});break
area=len(cells)*step*step
out={'schema':'r074-single-corner-Borel-full-cap-and-inner-cell-v1','n':n,'density':str(c),'unit_origin_cuts':[str(v)for v in cuts],'unit_vertex_count':len(vertices),'unit_vertices':vertices,'sharp_family_density':str(cap),'binding_unit_origin':[str(cu),str(cv)],'all_unit_vertices_passed':True,'receiver_mesh':N,'total_receiver_cells':N*N,'executed_receiver_event_tests':tests,'accepted_closed_cells':len(cells),'lower_area':str(area),'lower_area_decimal':float(area),'aux_corner_mass_over2':str(F(q*q,2*n*n)),'lower_minus_aux':str(area-F(q*q,2*n*n)),'main_threshold':str(F(1,8)+F(3,4*n**4)),'cell_certificate':cells,'scope':'Exact cap vertices plus strict whole-receiver-cell inner area certificate for density-c Borel nu=1_(P_SW minus J). General vertex-completeness and unit-to-all-squares proofs are separate. Hnu empty because density<=1. No finite atomic transfer.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.time()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k not in ['unit_vertices','cell_certificate','unit_origin_cuts']}))
