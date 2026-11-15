"""Read-only exact all-label contract check for an affine mass rectangle.
No geometry enumeration, no optimization, no imported envelope kernel.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,sys,time
start=time.process_time();path=Path(sys.argv[1]);target=Path(sys.argv[2]);out=Path(sys.argv[3]);data=json.loads(path.read_text());reference=json.loads(target.read_text())
base=data['base'];bx,by,lam=map(F,[base['Px'],base['Py'],base['lam']]);nx,ny=data['copy_factors'];Px,Py=nx*bx,ny*by
assert min(bx,by,lam)>0 and type(nx)is int and type(ny)is int and nx>0 and ny>0
A0,A1=map(F,data['a_bounds']);B0,B1=map(F,data['b_bounds']);assert A0<A1 and B0<B1
corners=data['corners'];assert len(corners)==4
assert all((F(c['Px']),F(c['Py']),F(c['lam']))==(Px,Py,lam) for c in corners)
n=len(base['pts']);assert all(len(c['pts'])==nx*ny*n for c in corners)
minimum=None;totals=[F(0)]*4;sum_u=sum_v=F(0);kappa=distance=F(0);records=[]
for idx,rows in enumerate(zip(*(c['pts'] for c in corners))):
 a=idx//(ny*n);b=(idx//n)%ny;i=idx%n
 x0,y0,w0=map(F,base['pts'][i]);assert w0>0;x0%=bx;y0%=by
 coords=[(F(z[0])%Px,F(z[1])%Py) for z in rows];assert len(set(coords))==1
 weights=list(F(z[2]) for z in rows);assert min(weights)>0 and weights[0]+weights[3]==weights[1]+weights[2]
 minimum=min(weights) if minimum is None else min(minimum,*weights)
 u=(weights[1]-weights[0])/(A1-A0);v=(weights[2]-weights[0])/(B1-B0);intercept=weights[0]-u*A0-v*B0;sum_u+=u;sum_v+=v
 for j,w in enumerate(weights):totals[j]+=w;kappa=max(kappa,w/w0)
 x,y=coords[0];zx,zy=x0+a*bx,y0+b*by;kx=((zx-x)/Px+F(1,2)).__floor__();ky=((zy-y)/Py+F(1,2)).__floor__();dx=x+kx*Px-zx;dy=y+ky*Py-zy
 distance=max(distance,abs(dx),abs(dy));assert max(abs(dx),abs(dy))<=F(data['delta'])
 records.append({'label':[a,b,i],'coordinate':list(map(str,coords[0])),'weights':list(map(str,weights)),'intercept':str(intercept),'u':str(u),'v':str(v),'lift':[kx,ky],'displacement':[str(dx),str(dy)]})
assert len(set(totals))==1 and sum_u==sum_v==0 and totals[0]>0
base_mass=sum(F(z[2]) for z in base['pts']);delta=F(data['delta']);R=F(5)
tail=kappa*base_mass/(bx*by)*(1+(bx+2*delta)/(2*R))*(1+(by+2*delta)/(2*R));assert tail<lam
# The exact R20 source is at the saved point (-3/20,33/4); match every label.
param_a,param_b=F(-3,20),F(33,4);assert A0<=param_a<=A1 and B0<=param_b<=B1
assert (F(reference['Px']),F(reference['Py']),F(reference['lam']))==(Px,Py,lam) and len(reference['pts'])==len(records)
for row,record in zip(reference['pts'],records):
 x,y,w=map(F,row);assert [str(x%Px),str(y%Py)]==record['coordinate']
 assert w==F(record['intercept'])+param_a*F(record['u'])+param_b*F(record['v'])
result={'scope':'exact R20 all-label family contract only; no halo geometry',
 'python':platform.python_version(),'cpu_seconds':time.process_time()-start,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'target_input_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'bounds':{'a':list(map(str,[A0,A1])),'b':list(map(str,[B0,B1]))},'labels':len(records),'copy_factors':[nx,ny],'all_corners_strictly_positive':True,'minimum_corner_mass':str(minimum),'actual_constant_mass':str(totals[0]),'sum_u':str(sum_u),'sum_v':str(sum_v),'uniform_kappa':str(kappa),'maximum_lift_distance':str(distance),'delta':str(delta),'R':str(R),'uniform_tail':str(tail),'margin':str(lam-tail),'reference_point':list(map(str,[param_a,param_b])),'R20_point_matches_every_label':True,'all_label_records':records}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='all_label_records'},indent=2))
