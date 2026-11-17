"""Independent Gaussian-erfc evaluation; supporting check of rational certificate."""
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp,json,os,sys,platform,hashlib
mp.mp.dps=90
base=Path(__file__).parent;certpath=base/'round013-core-results.json';cert=json.loads(certpath.read_text());rows=[]
def number(s):
 q=Q(s);return mp.mpf(q.numerator)/q.denominator
for row in cert['certificates']:
 src=base/row['input_file'];d=json.loads(src.read_text(),parse_float=str);knots=list(map(number,d['parameters']['knots']));vals=list(map(number,d['best']['values']));lam=number(row['lambda_value']);total=mp.mpf(0)
 for l,r,fl,fr in zip(knots,knots[1:],vals,vals[1:]):
  for aa,bb,fa,fb in [(l,r,fl,fr),(-r,-l,fr,fl)]:
   m=(aa+bb)/2;radius=(bb-aa)/2;slope=(fb-fa)/(bb-aa);fm=(fa+fb)/2
   q0=-m*m/2+lam*m-lam*lam*fm/2;a=lam-m-lam*lam*slope/2
   u=(-radius-a)/mp.sqrt(2);v=(radius-a)/mp.sqrt(2)
   if u>0:diff=mp.erfc(u)-mp.erfc(v)
   elif v<0:diff=mp.erfc(-v)-mp.erfc(-u)
   else:diff=mp.erf(v)-mp.erf(u)
   total+=mp.exp(q0+a*a/2)*diff/2
 lo=number(row['core_lower_dyadic']);hi=number(row['core_upper_dyadic']);assert lo<=total<=hi
 rows.append(dict(input=row['input_file'],core_erfc=str(total),inside_rational_interval=True,interval_width=str(hi-lo)))
out=dict(rows=rows,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,dps=mp.mp.dps,pid=os.getpid()),certificate_sha256=hashlib.sha256(certpath.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='90-digit independent numerical crosscheck; rigorous status belongs to rational certificate')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
