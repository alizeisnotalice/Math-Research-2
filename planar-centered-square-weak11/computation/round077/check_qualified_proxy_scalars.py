#!/usr/bin/env python3
"""Exact scalar and finite-cell checks; no full halo integration or cap enumeration."""
from fractions import Fraction as F
from pathlib import Path
import json,sys,os,time,platform,hashlib
start=time.perf_counter();ip=Path(sys.argv[1]);cfg=json.loads(ip.read_text())
n=cfg['n'];e=F(1,n*n);h=F(n+1,n*n);w=1-e;m=F((n//2)**2,n*n)
a,d,b1,b2=(F(cfg[k])for k in ['a','d','b1','b2']);u2=(m+b1)/4;v2=(m+b2)/4
assert n==20 and m==b1==F(1,4) and b2==F(1,2) and d==a+F(1,2)
assert u2==F(1,8) and v2==F(3,16) and m+b1+b2==1
# Every comparison involving a radical is squared only after checking positivity.
assert F(1,4)<a-e and (a-e)**2<u2<a*a<v2<F(1,4)
assert d*d<4*v2 and d-F(1,2)>0 and (d-F(1,2))**2>u2
assert (a-F(1,4))**2<u2 and F(1,4)>e
# Expand the two disjoint true-residual area contributions as coeffs of 1,u,v.
# (u-e)^2-(u-a+1/4)(1/4-e)+(1/2-e)(1/2-a)-(v-e)(2v-d)
constant=u2+e*e+(a-F(1,4))*(F(1,4)-e)+(F(1,2)-e)*(F(1,2)-a)-2*v2-d*e
coef_u=-2*e-(F(1,4)-e);coef_v=d+2*e
assert constant==-F(24641,160000) and coef_u/4==-F(101,1600) and coef_v/4==F(43,200)
r2=F(cfg['sqrt2_upper']);r3=F(cfg['sqrt3_lower'])
assert r2>0 and r2*r2>2 and r3>0 and r3*r3<3
lower=constant+coef_u*r2/4+coef_v*r3/4
assert lower==F(206483,1600000) and lower-m/2==F(6483,1600000)>0
lo,hi=map(F,cfg['receiver_cell_positive']);r=F(cfg['proxy_witness_radius'])
assert e<lo<hi<F(1,2) and hi<r<F(1,2)
# Each collapsed source is strictly captured on the whole closed cell.
collapsed=[(F(0),F(0)),(-a,F(0)),(F(0),-d)]
maxdistance=max(abs(x-z)for q in collapsed for z in q for x in [-hi,-lo])
assert maxdistance==hi<r and 1-4*r*r==F(199,10000)>0
assert F(1,2)-lo<h and e+b1+b2<4*lo*lo
assert b1+b2<4*lo*lo and lo>F(1,4)
# Unit-window captures before and after replacing local old mass by m at zero.
Q=(-a,1-a,-d,1-d);old=[(i*h,j*h)for i in range(n)for j in range(n)]
inQ=lambda z: Q[0]<=z[0]<=Q[1] and Q[2]<=z[1]<=Q[3]
inP=lambda z: -1<z[0]<F(1,2)-e and -1<z[1]<F(1,2)-e
oldcount=sum(map(inQ,old));remaining=sum(inQ(z)and not inP(z)for z in old)
assert sum(map(inP,old))==100 and oldcount==39 and remaining==9
assert oldcount*e+b1+b2==F(339,400) and remaining*e+m+b1+b2==F(409,400)
# Explicit strict eligible edges connect the entire actual support.
assert 2*e>h*h and e+b1>a*a and b1+b2>d*d
nonterminal=(w+d)**2-F(35,33)*(1+b1+b2)
assert nonterminal==F(8319673,5280000)>0
out={'passed':True,'proxy_exact':'43*sqrt(3)/200-101*sqrt(2)/1600-24641/160000','proxy_strict_rational_lower':str(lower),'strict_half_local_mass_excess_lower':str(lower-m/2),'proxy_only_cell_area':str((hi-lo)**2),'cell_common_open_radius':str(r),'cell_proxy_mass_margin':str(1-4*r*r),'cell_actual_small_radius_mass_upper':str(e+b1+b2),'cell_threshold_lower':str(4*lo*lo),'old_unit_capture_count':oldcount,'postcollapse_retained_old_count':remaining,'original_unit_mass':str(oldcount*e+b1+b2),'postcollapse_unit_mass':str(remaining*e+m+b1+b2),'nonterminal_gap':str(nonterminal),'scope':'Exact scalar identities, actual 400-point unit-window count, explicit receiver-cell separation, and connecting-edge inequalities. No full global cap enumeration or halo area integration; those belong to the separate production certificate.','runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.perf_counter()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest()}}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
