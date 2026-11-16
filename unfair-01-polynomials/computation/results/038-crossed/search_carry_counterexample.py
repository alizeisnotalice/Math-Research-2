import json,os,time
from pathlib import Path
assert os.path.ismount('/Volumes/data')
def mask(weights):
 s=1
 for w in weights:
  if s&(s<<w):return None
  s|=s<<w
 return s
best=None
for b in range(3,22,2):
 s=mask([8,12,18,27,12*b,18*b])
 if s is None:continue
 for L in range(1,1000,3):
  if s&(s<<(2*L)):continue
  z=s|(s<<(2*L))
  if z&(z<<(3*L)):continue
  z|=z<<(3*L)
  for M in range(2,1000,3):
   total=6*b+13+L+M
   if best and total>=best['H_degree']:continue
   if z&(z<<(2*M)):continue
   t=z|(z<<(2*M))
   if t&(t<<(3*M)):continue
   t|=t<<(3*M)
   best={'b':b,'L':L,'M':M,'H_degree':total,'C_degree':t.bit_length()-1,'C_terms':t.bit_count()}
print(best)
Path(__file__).with_suffix('.json').write_text(json.dumps(best,indent=2)+'\n')
