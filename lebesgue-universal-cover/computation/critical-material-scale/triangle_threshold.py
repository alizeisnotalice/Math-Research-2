from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json,sys
from flint import arb,ctx
ctx.prec=256;comp=Path(sys.argv[1]).resolve();sys.path.insert(0,str(comp))
from certify_polygon_083731 import verified_table,Q,ar
from certify_tomar_six_high import SourceTree
tree=SourceTree(Path(sys.argv[2]));source=json.loads((comp/'tomar-six-high-input.json').read_text());s0=F('.99999900555134042347742929191955');s1=F('.999999005551341');s=ar(s0);pi=arb.pi();r=s/(2*(pi/6).cos());out=[];a=arb(81827130)/Q;c=(1+1/arb(3).sqrt())*a;critical=(c+(c*c+arb(16)*a*a/3).sqrt())*arb(3)/4;assert critical>ar(s0) and critical<ar(s1)
for idx,level in ((3,32),(4,64)):
 path=Path(sys.argv[idx]);counts,points,_=verified_table(path,tree,level,source['tables'][str(level)]);p=points[0,counts[0]:counts[0]+counts[1]];bad=[];minok=None
 for j,(x,y) in enumerate(p):
  X,Y=ar(F(int(x),Q)),ar(F(int(y),Q))
  for k in range(3):
   sn,co=(2*pi*k/3).sin_cos();g=s*s-(X-r*co)**2-(Y-r*sn)**2
   if not g>0:bad.append([j,k,int(x),int(y),str(g)])
   elif minok is None or g.lower()<minok:minok=g.lower()
 expect=[[12,2,81827130,81827130],[84,1,81827130,-81827130]] if level==32 else [[24,2,81827130,81827130],[168,1,81827130,-81827130]]
 assert [row[:4] for row in bad]==expect and minok>0
 out.append({'table':level,'not_positive':bad,'min_other_squared_margin':str(minok),'table_sha256':sha256(path.read_bytes()).hexdigest()})
Path(sys.argv[5]).write_text(json.dumps({'lower_s':str(s0),'upper_s':str(s1),'critical':str(critical),'results':out},indent=2)+'\n')
