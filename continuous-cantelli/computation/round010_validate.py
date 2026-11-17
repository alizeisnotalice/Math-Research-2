"""Independent high-precision integration of frozen floating coefficients."""
from pathlib import Path
import mpmath as mp,json,sys,os,platform,hashlib
mp.mp.dps=60
src=Path(__file__).with_name('round010-finite-results.json');d=json.loads(src.read_text());row=d['rows'][-1]
a=list(map(lambda x:mp.mpf(str(x)),row['coefficients']));centers=list(map(lambda x:mp.mpf(str(x)),d['parameters']['centers']));width=mp.mpf(str(d['parameters']['width']))
def integrate(kind,lam=0):
 result=mp.mpf(0)
 for center,coef in zip(centers,a):
  def fun(y):
   if abs(y)>=1:return mp.mpf(0)
   x=center+width*y;psi=mp.exp(1-1/(1-y*y));b=coef*psi
   weight=2*width*mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
   if kind=='mean':v=b
   elif kind=='fourth':v=b*b-2*(x*x-1)*b
   else:v=mp.exp(-lam*lam/2)*mp.cosh(lam*x)*mp.expm1(-lam*lam*b/2)
   return weight*v
  result+=mp.quad(fun,[-1,mp.mpf('-.5'),0,mp.mpf('.5'),1])
 return str(result)
out=dict(mean=integrate('mean'),fourth_condition=integrate('fourth'),F_minus_one={str(l):integrate('F',mp.mpf(l)) for l in [1,2,3,4]},runtime=dict(python=platform.python_version(),mpmath=mp.__version__,dps=mp.mp.dps,pid=os.getpid()),input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='60-digit quadrature of decimalized floating coefficients; not interval certification or exact zero')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
