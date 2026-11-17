"""Rational certificates for uniform compact-tail detection and tail leakage."""
from fractions import Fraction as Q
from math import factorial,ceil
from pathlib import Path
import json,sys,os,platform,hashlib
import mpmath as mp
mp.mp.dps=70
rows=[]
for R,L,d in [(4,11,Q(2,5)),(8,11,Q(2,5)),(12,11,Q(2,5)),(8,1,Q(2,5)),(8,100,Q(2,5)),(8,11,Q(1,5))]:
 for lam in range(R+1,2001):
  r=min(Q(1),d/(2*L),Q(1,L*lam*lam));h=d-L*r;m0=2*r/5
  E=h*lam*lam/2-(R+r)*lam+Q(R*R,2)
  if h*lam*lam/2<1 or E<=0:continue
  lower=sum(E**j/factorial(j) for j in range(33))
  if lower<4/m0:continue
  A=(1-h)*lam*lam/2+r*lam;n=ceil(A);pref=m0/4
  assert 0<r<=1 and h>0 and h*lam*lam/2>=1
  assert E==h*lam*lam/2-(R+r)*lam+Q(R*R,2)
  assert lower>=4/m0 and n>=A
  logeta=mp.log10(mp.mpf(pref.numerator)/pref.denominator)-n*mp.log10(3)
  rows.append(dict(R=R,L=L,delta=str(d),lambda_value=lam,r=str(r),height=str(h),mass_lower=str(m0),ratio_exponent=str(E),exponential_series_degree=32,series_lower=str(lower),required_series_lower=str(4/m0),A=str(A),ceil_A=n,eta_factor=str(pref),eta_formula=f'({pref}) * 3^(-{n})',tail_integral_factor=str(2*pref/(lam*lam)),log10_eta_display=str(logeta)))
  break
 else:raise RuntimeError('no certificate within search range')
applications=[]
for filename in ['round004-L8-T12-long.json','round004-L12-T12.json']:
 src=Path(__file__).with_name(filename);dd=json.loads(src.read_text(),parse_float=str)
 knots=list(map(Q,dd['parameters']['knots']));values=list(map(Q,dd['best']['values']))
 assert knots[0]==0 and values[0]==Q(3,5) and values[-1]==1
 slope=max(abs((values[j+1]-values[j])/(knots[j+1]-knots[j])) for j in range(len(knots)-1))
 assert slope<=11 and min(values)>0
 cert=next(row for row in rows if row['R']==knots[-1] and row['L']==11 and row['delta']=='2/5')
 applications.append(dict(input_file=filename,input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),slope_exact=str(slope),anchor=str(values[0]),tail_convention='even linear interpolation of knots; f=1 outside last knot, as defined in round4',lambda_value=cert['lambda_value'],eta_formula=cert['eta_formula']))
out=dict(certificates=rows,applications=applications,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,dps=mp.mp.dps,pid=os.getpid()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='all certificate inequalities are exact Fraction arithmetic; log10 displays only are non-directed high precision')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
for r in rows: print({k:r[k] for k in ['R','L','delta','lambda_value','eta_formula','log10_eta_display']})
