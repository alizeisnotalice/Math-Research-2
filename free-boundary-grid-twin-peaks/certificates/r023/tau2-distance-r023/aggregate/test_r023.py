#!/usr/bin/env python3
import json,hashlib,pathlib,subprocess,sys,tempfile
from fractions import Fraction
ROOT=pathlib.Path(__file__).parent; A=ROOT/'aggregate_tau2_distance.py'; RAW=ROOT/'n3-p17-r2-distance.jsonl'; SRC=pathlib.Path('/Volumes/data/CodexTemp/T2/a022-tau2-distance-terminal-r023/tau2_distance_terminal_r023.swift')
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def run(raw,side,ok,eps='1/3'):
 c=subprocess.run([sys.executable,str(A),str(raw),'--sidecar',str(side),'--source',str(SRC),'--epsilon',eps],text=True,capture_output=True)
 if (c.returncode==0)!=ok:raise AssertionError((ok,c.stderr))
 return json.loads(c.stdout) if ok else None
def main():
 def inds(n,e):
  q=Fraction(e);return [d for d in range(2,2*n-1) if d*q.denominator>n*q.numerator]
 assert inds(7,'2/7')[0]==3 and inds(7,'3/7')[0]==4
 def indsF(n,q):return [d for d in range(2,2*n-1) if d*q.denominator>n*q.numerator]
 assert indsF(7,Fraction(2,7)-Fraction(1,10**90))[0]==2 and indsF(7,Fraction(2,7)+Fraction(1,10**90))[0]==3
 # Independent hand oracle specified in the review: Z=(1,3), L=(4,4),
 # L2=(10,8), hist=(1,3,0),(2,0,2), all prior to invoking the program.
 z=[Fraction(1),Fraction(3)];l=[Fraction(4),Fraction(4)];l2=[Fraction(10),Fraction(8)];h=[[Fraction(1),Fraction(3),Fraction(0)],[Fraction(2),Fraction(0),Fraction(2)]];S=sum(zz*ll for zz,ll in zip(z,l));p=[sum(z[r]*h[r][d] for r in range(2))/S for d in range(3)];e=sum(z[r]*(l2[r]+3*l[r]) for r in range(2))/(2*S)
 def se2(x,fac=1):return Fraction(2,1)*sum(t*t for t in x)/S/S/fac/fac
 assert p==[Fraction(7,16),Fraction(3,16),Fraction(3,8)] and e==Fraction(41,16)
 assert se2([z[r]*(h[r][0]-p[0]*l[r]) for r in range(2)])==Fraction(9,1024)
 assert se2([z[r]*(h[r][1]-p[1]*l[r]) for r in range(2)])==Fraction(81,1024)
 assert se2([z[r]*(h[r][2]-p[2]*l[r]) for r in range(2)])==Fraction(9,256)
 assert se2([z[r]*(l2[r]+3*l[r]-2*e*l[r]) for r in range(2)],2)==Fraction(9,1024)
 rows=[json.loads(x) for x in RAW.read_text().splitlines()]
 with tempfile.TemporaryDirectory() as td:
  td=pathlib.Path(td)
  def put(name,x):
   r=td/(name+'.jsonl');r.write_text('\n'.join(json.dumps(y) for y in x)+'\n');m=x[0];s=td/(name+'.json');json.dump({'raw_sha256':sha(r),'source_sha256':sha(SRC),'expected_metadata':{k:m[k] for k in ('n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps')}},open(s,'w'));return r,s
  r,s=put('good',rows);o=run(r,s,True);assert o['strict_tail_d_gt_epsilon_n']=='1.0000000000000000000000000000000000000000000000000000000000000000000000000000000'
  # strict d>2: tail equals 1-p(d=2); strict d>3 equals p(d=4).
  a=run(r,s,True,'2/3');b=run(r,s,True,'1');p=o['distance_L_weighted_probabilities'];assert abs(float(a['strict_tail_d_gt_epsilon_n'])-(1-float(p[2])))<1e-12 and abs(float(b['strict_tail_d_gt_epsilon_n'])-float(p[4]))<1e-12
  x=json.loads(json.dumps(rows));x[1]['tau2_distance_L_hist_groups'][0][2]+=1;r,s=put('hist',x);run(r,s,False)
  x=json.loads(json.dumps(rows));x[1]['tau2_distance_L2_groups'][0]-=1;r,s=put('l2',x);run(r,s,False)
  x=json.loads(json.dumps(rows));z=x[1];z.update(Z=0,logZ=-1000,tau2_lift_complete_K2_Z=0,tau2_lift_complete_K2_D=0,tau2_lift_complete_K2_F=0,tau2_lift_complete_K2_G=0);r,s=put('underflow',x);run(r,s,True)
 print('r023_tests=PASS: unequal actual R2 ratio/delta path, hist, L2, strict epsilon boundaries, underflow')
if __name__=='__main__':main()
