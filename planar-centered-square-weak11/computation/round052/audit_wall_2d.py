#!/usr/bin/env python3
"""Rational enclosure audit: all active bbox intervals near a candidate wall, z fixed."""
import json, hashlib, sys, gzip
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
p=Path(__file__).resolve().parent
def slab_path(n):
 q=p/n/'slab.json'; gz=q.with_suffix('.json.gz'); return gz if gz.exists() else q
def slab_read(n):
 q=slab_path(n); return json.loads(gzip.open(q,'rt').read() if q.suffix=='.gz' else q.read_text())
D=[slab_read(n) for n in ('one_axis_t0','one_axis_t1')]
a0=json.loads((p/'slice/corner_t0_z0.json').read_text()); a1=json.loads((p/'slice/corner_t1_z0.json').read_text()); v=json.loads((p/'slice/velocity_mass.json').read_text())
assert len(a0['pts'])==len(a1['pts'])==len(v) and all(a0[k]==a1[k] for k in a0 if k!='pts')
assert all(a[:2]==b[:2] and F(b[2])-F(a[2])==F(w) and F(a[2])>0 and F(b[2])>0 for a,b,w in zip(a0['pts'],a1['pts'],v))
assert all(D[0][k]==D[1][k] for k in ('x_classes','y_classes','domain','s0','periodic','source_order_sha256'))
assert len(D[0]['source_catalog'])==len(D[1]['source_catalog'])
assert all({k:v for k,v in a.items() if k!='w'}=={k:v for k,v in b.items() if k!='w'} for a,b in zip(D[0]['source_catalog'],D[1]['source_catalog']))
W=json.loads((p/'wall_witness.json').read_text()); T=[F(W['t_lower'])-F(1,10**9),F(W['t_upper'])+F(1,10**9)]; X=[F(W['x_lower'])-F(1,10**10),F(W['x_upper'])+F(1,10**10)]; Z=F(1,10**11); assert all(F(d['H'])>=Z and F(d['delta_empty'])>=Z for d in D); Y=F(W['common_y']); YI=[Y-F(1,10**7),Y+F(1,10**7)]; N=1<<128
m=[{tuple(r['bbox']):b for r,b in zip(d['box_records'],d['boxes'])} for d in D]; assert m[0].keys()==m[1].keys()
assert all(F(d['tail']['radius'])==10 for d in D)
for k in m[0]:
 a,b=m[0][k],m[1][k]
 assert all(e[:3]==f[:3] and F(e[3])==F(a[0][3]) and F(f[3])==F(b[0][3]) for e,f in zip(a,b))
def sqrt_bounds(q):
 k=isqrt(q.numerator*N*N//q.denominator); return F(k,N),F(k+(k*k*q.denominator!=q.numerator*N*N),N)
def ends(k):
 a,b=m[0][k],m[1][k];q0,q1=F(a[0][3]),F(b[0][3]);q=[q0+(q1-q0)*t for t in T];rl=sqrt_bounds(min(q))[0];rh=sqrt_bounds(max(q))[1]
 return [(F(e[0])-rh-abs(F(e[1]))*Z,F(e[0])-rl+abs(F(e[1]))*Z) if e[2]==-1 else (F(e[0])+rl-abs(F(e[1]))*Z,F(e[0])+rh+abs(F(e[1]))*Z) if e[2]==1 else (F(e[0])-abs(F(e[1]))*Z,F(e[0])+abs(F(e[1]))*Z) for e in a]
possible=[]; margins=[]
for k in m[0]:
 e=ends(k); sep=[e[0][0]-X[1],X[0]-e[1][1],e[2][0]-YI[1],YI[0]-e[3][1]]; v=max(sep)
 if v>0:margins.append(v)
 else:possible.append(k)
# Independently verify the stored radical crossing interval from endpoint input data.
a0b=m[0][(86,88,29,29)][0]; b0b=m[0][(85,88,29,29)][1]; b1b=m[1][(85,88,29,29)][1]
assert F(a0b[1])==0 and F(b0b[1])==1 and a0b==m[1][(86,88,29,29)][0]
ra0,ra1=sqrt_bounds(F(a0b[3])); c=F(a0b[0])-F(b0b[0]); assert c>ra1
q0,q1=F(b0b[3]),F(b1b[3]); assert q1<q0
root=sorted(((c-ra0)**2-q0)/(q1-q0) for ra0 in (ra0,ra1))
assert F(W['t_lower'])<=root[0] and root[1]<=F(W['t_upper'])
expected=[(85,88,29,29),(86,88,29,29),(87,88,29,29),(88,88,29,29)]; assert possible==expected
a,b,c,d=expected
for j in (0,1):
 assert all(m[j][k][i]==m[j][b][i] for k in (c,d) for i in (0,2,3))
for k in expected:
 e=ends(k); assert e[2][1]<YI[0] and e[3][0]>YI[1]
assert ends(a)[0][1]<X[0] and all(ends(k)[1][0]>X[1] for k in (b,c,d))
assert all(X[0]<ends(k)[i][0] and ends(k)[i][1]<X[1] for k,i in ((a,1),(b,0)))
print('potential boxes',len(possible),possible,'excluded',len(margins),'minimum exclusion',float(min(margins)))
out={'scope':'On the specified physical position interval around physical z=P8/1000000 (P8 is half the input supercell period Px), endpoint slab omitted-box radii and affine masses ensure the same bbox list contains every potentially active box throughout t in [0,1]. Rational interval exclusions hold on the stated t,x,y,position neighborhood.','position_s0':D[0]['s0'],'position_half_width':str(Z),'t_interval':list(map(str,T)),'x_interval':list(map(str,X)),'y_interval':list(map(str,YI)),'possible_boxes':[list(k) for k in possible],'excluded_count':len(margins),'minimum_exclusion_margin':str(min(margins)),'input_sha256':{n:hashlib.sha256((p/'slice'/n).read_bytes()).hexdigest() for n in ('corner_t0_z0.json','corner_t1_z0.json','velocity_mass.json')},'slab_sha256':[hashlib.sha256(gzip.open(slab_path(n),'rb').read() if slab_path(n).suffix=='.gz' else slab_path(n).read_bytes()).hexdigest() for n in ('one_axis_t0','one_axis_t1')],'witness_sha256':hashlib.sha256((p/'wall_witness.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version}
f=p/'wall_local_2d.json';f.write_text(json.dumps(out,indent=2)+'\n');print('sha256',hashlib.sha256(f.read_bytes()).hexdigest())
