"""Exact arithmetic audit of degree-seven envelopes and the finite support case split."""
import argparse,itertools,json,math,os,sys
from fractions import Fraction as F
from pathlib import Path

def main():
 if not __debug__:raise RuntimeError('no -O')
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists()
 c=[F(1)]
 for r in [F(2)]*3+[F(1,2)]*3+[F(1)]:
  b=[F(0)]*(len(c)+1)
  for i,v in enumerate(c):b[i]+=v;b[i+1]+=r*v
  c=b
 assert c==list(map(F,[1]))+[F(17,2),F(117,4),F(419,8),F(419,8),F(117,4),F(17,2),F(1)]
 powers=[3*(F(2)**j+F(1,2)**j)+1 for j in range(1,8)]
 assert powers==[F(17,2),F(55,4),F(203,8),F(787,16),F(3107,32),F(12355,64),F(49283,128)]
 allowed=[]
 for n in range(1,7):
  for S in itertools.combinations(range(1,7),n):
   r,t=S[0],7-S[-1]
   forbidden=r==t or (r==1 and t in (2,3,4)) or (t==1 and r in (2,3,4)) or (r,t) in ((2,4),(4,2))
   if not forbidden:allowed.append(list(S))
 wanted=[[1],[2],[3],[4],[5],[6],[1,2],[2,4],[3,5],[5,6],[2,3,4],[3,4,5]]
 assert sorted(allowed)==sorted(wanted)
 fields=[D for D in range(2,190) if all(D%(i*i) for i in range(2,math.isqrt(D)+1)) and (D%4==1 or D<=47)];assert len(fields)==58
 report={'status':'passed','coefficient_envelope':list(map(str,c)),'newton_bounds':list(map(str,powers)),'surviving_nonrational_supports':allowed,'fields':fields,'scope':'Exact constants and finite case split only; analytic envelope and endpoint assumptions require the accompanying proof.','python':sys.version}
 t=a.out.with_suffix('.partial');t.write_text(json.dumps(report,indent=2)+'\n');os.replace(t,a.out);print(json.dumps(report))
if __name__=='__main__':main()
