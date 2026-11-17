"""Independent mpmath integration of one interval witness. Not interval arithmetic."""
from pathlib import Path
import json,sys,os,platform,hashlib
import mpmath as mp
mp.mp.dps=60
src=Path(sys.argv[1]);d=json.loads(src.read_text(),parse_float=str)
k=list(map(mp.mpf,d['parameters']['knots']));v=list(map(mp.mpf,d['best']['values']));L=k[-1]
t=mp.mpf('1.7');left=mp.mpf('0.475');right=mp.mpf('0.525')
phi=lambda y:mp.exp(-y*y/2)/mp.sqrt(2*mp.pi)
Phi=lambda y:mp.erfc(-y/mp.sqrt(2))/2
def diff(lo,hi):
 return (mp.erfc(lo/mp.sqrt(2))-mp.erfc(hi/mp.sqrt(2)))/2 if lo>0 else Phi(hi)-Phi(lo)
def kernel(y,var):
 sd=mp.sqrt(var)
 return phi(y)*(diff((left-y)/sd,(right-y)/sd)+diff((left+y)/sd,(right+y)/sd))
source=mp.mpf(0);segments=0
for i in range(len(k)-1):
 a,b=k[i:i+2];fa,fb=v[i:i+2];slope=(fb-fa)/(b-a)
 if fa>=t and fb>=t:continue
 if fa<t and fb<t:lo,hi=a,b
 elif slope>0:lo,hi=a,a+(t-fa)/slope
 else:lo,hi=a+(t-fa)/slope,b
 if hi<=lo:continue
 points=sorted(set([lo,hi]+[q for q in [left,right,(lo+hi)/2] if lo<q<hi]))
 def integrand(y):
  var=t-fa-slope*(y-a)
  return mp.mpf(0) if var<=0 else kernel(y,var)
 source+=mp.quad(integrand,points);segments+=1
source+=mp.quad(lambda y:kernel(y,t-1),[L,L+2,mp.inf])
free=diff(left/mp.sqrt(t),right/mp.sqrt(t));mass=free-source
out={'input':src.name,'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
 'dps':mp.mp.dps,'t':str(t),'interval':[str(left),str(right)],'mass':mp.nstr(mass,55),
 'free_mass':mp.nstr(free,55),'source_mass':mp.nstr(source,55),'active_segments':segments,
 'runtime':{'python':platform.python_version(),'mpmath':mp.__version__,'pid':os.getpid()},
 'scope':'high-precision floating evaluation, not directed interval certificate',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
