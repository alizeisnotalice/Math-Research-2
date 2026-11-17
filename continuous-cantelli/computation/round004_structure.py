"""Exact structural audit of decimal knot values; not a transform certificate."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,math,hashlib,os,platform
rows=[]
for path in sys.argv[2:]:
 p=Path(path);d=json.loads(p.read_text(),parse_float=str)
 k=list(map(F,d['parameters']['knots']));v=list(map(F,d['best']['values']));L=k[-1]
 assert k[0]==0 and all(k[i+1]-k[i]==F(1,4) for i in range(len(k)-1))
 assert v[0]==v[1]==F(3,5) and v[-1]==1
 assert min(v)>=F(1,20) and max(v)<=F(14,5)
 slope=max(abs((v[i+1]-v[i])/(k[i+1]-k[i])) for i in range(len(k)-1))
 assert slope<=11
 assert L>=1
 lam=5*(L+1)
 margin=F(15,4)*(L+1)
 # gamma([-1/4,1/4]) > (1/2)*(1/3)*(31/32)>4/25, using pi<4.
 mass_lower=F(31,192)
 assert mass_lower>F(4,25)
 # For L>=1, (2/25)*exp(margin)-1 >= (2/25)*(1+15/2+(15/2)^2/2)-1.
 bracket_lower=F(2,25)*(1+F(15,2)+F(15,2)**2/2)-1
 assert bracket_lower==F(193,100)>0
 log_bound=float(-lam*lam/2+L*lam)+math.log(2/25)+float(margin)+math.log1p(-12.5*math.exp(-float(margin)))
 rows.append(dict(input_file=p.name,input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),L=str(L),
                  minimum_f=str(min(v)),maximum_f=str(max(v)),exact_maximum_slope=str(slope),
                  uniform_lipschitz_bound='11',plateau_value='3/5',plateau_radius='1/4',
                  tail_value='1',witness_lambda=str(lam),positive_bracket_lower=str(bracket_lower),
                  lower_bound_formula='exp(-lambda^2/2+L*lambda)*((2/25)*exp(15*(L+1)/4)-1)',
                  lower_bound_log10_approx=log_bound/math.log(10),
                  strict_exclusion='finite tail and min f>0 force any Gaussian output variance to H=3; witness F(lambda)>1 excludes it'))
out=dict(certificates=rows,arithmetic='exact Fraction for structural assertions; logarithmic display is float only',
         runtime=dict(python=platform.python_version(),pid=os.getpid()),
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k not in ['minimum_f','maximum_f','exact_maximum_slope','input_sha256']} for r in rows],indent=2))
