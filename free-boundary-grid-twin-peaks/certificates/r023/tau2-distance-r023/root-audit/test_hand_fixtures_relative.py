"""Root integration test: known rational fixture is submitted to the actual CLI."""
import copy,hashlib,json,math,pathlib,subprocess,sys,tempfile
from decimal import Decimal as D
from fractions import Fraction as F
ROOT=pathlib.Path(__file__).parent
HERE=pathlib.Path(__file__).resolve().parent
AGG=HERE.parent/'aggregate/aggregate_tau2_distance.py'
SRC=HERE.parent/'terminal/tau2_distance_terminal_r023.swift'
TEMPLATE=HERE.parent/'aggregate/n3-p17-r2-distance.jsonl'
base=[json.loads(s) for s in TEMPLATE.read_text().splitlines()]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fixture(n):
 m=copy.deepcopy(base[0]);m.update(n=n,particles=2,replicates=2,tau_lo=2,tau_hi=2,domain='0',move_steps=0)
 out=[m]
 for r,(Z,Ls,ds) in enumerate([(1,[1,3],[2,3]),(3,[2,2],[2,4])]):
  rec=copy.deepcopy(base[1]);rec.update(replicate=r,n=n,tau=2,dead=False,seed_hex=f'{0xA022000000000000+n*1000000+20000+r:016x}',Z=Z,logZ=math.log(Z))
  de=[];ff=[];gg=[];hist=[0]*(2*n-1)
  for L,d in zip(Ls,ds):
   x=(0,0);y={2:(1,1),3:(2,1),4:(2,2)}[d];t=[-(n-1),-(n-1)];u=[2*y[0]-(n-1),2*y[1]-(n-1)]
   de.append(y[0]**2+y[1]**2);ff.append(sum(v*v for v in t+u));gg.append(sum(t[i]*u[i] for i in range(2)));hist[d]+=L
  w=sum(Ls);ld=sum(L*v for L,v in zip(Ls,de));lf=sum(L*v for L,v in zip(Ls,ff));lg=sum(L*v for L,v in zip(Ls,gg));l2=sum(L*L for L in Ls)
  assert lf-2*lg==4*ld
  bd=sum(de)/2;bf=sum(ff)/2;bg=sum(gg)/2
  rec.update(U=Z*bd/(n*n),U_over_Z=bd/(n*n),scaled_U=bd/(n*n),D=Z*bd,F=Z*bf,G=Z*bg,F_over_Z=bf,G_over_Z=bg,V=bf/8,C=bg/4,F_minus_2G_minus_4D=0)
  for key,v in [('L',w),('L_distance2',ld),('L_F0',lf),('L_G0',lg)]:rec['tau2_lift_'+key+'_groups']=[v]
  rec['tau2_lift_group_LF_minus_2LG_minus_4LD']=[0]
  for key,v in [('L',w),('L_distance2',ld),('L_F0',lf),('L_G0',lg)]:rec['tau2_lift_scaled_mean_'+key]=v/2
  for key,v in [('Z',w),('D',ld),('F',lf),('G',lg)]:rec['tau2_lift_complete_K2_'+key]=Z*v/2
  rec.update(tau2_distance_L_hist_groups=[hist],tau2_distance_L2_groups=[l2],tau2_distance_scaled_mean_L_hist=[v/2 for v in hist],tau2_distance_scaled_mean_L2=l2/2)
  out.append(rec)
 out.append({'kind':'run_end','wall_seconds':0})
 return out

def near(actual,expected):
 target=D(expected.numerator)/D(expected.denominator)
 assert abs(D(actual)-target)<D('1e-12'),(actual,str(expected))
passed=[]
with tempfile.TemporaryDirectory(dir=HERE.parent) as td:
 td=pathlib.Path(td)
 def run(rows,eps,ok=True):
  raw=td/'fixture.jsonl';raw.write_text('\n'.join(json.dumps(x) for x in rows)+'\n');side=td/'side.json';m=rows[0]
  side.write_text(json.dumps({'raw_sha256':sha(raw),'source_sha256':sha(SRC),'expected_metadata':{k:m[k] for k in ('n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps')}}))
  c=subprocess.run([sys.executable,str(AGG),str(raw),'--sidecar',str(side),'--source',str(SRC),'--epsilon',str(eps)],text=True,capture_output=True)
  assert (c.returncode==0)==ok,c.stderr
  return json.loads(c.stdout) if ok else None
 for n in [3,7]:
  rows=fixture(n);o=run(rows,F(2,n))
  for d,p,v in [(2,F(7,16),F(9,1024)),(3,F(3,16),F(81,1024)),(4,F(3,8),F(9,256))]:
   near(o['distance_L_weighted_probabilities'][d],p);near(D(o['distance_deltaSE'][d])**2,v)
  near(o['E2_tau'],F(41,16));near(D(o['deltaSE_E2_tau'])**2,F(9,1024));near(o['strict_tail_d_gt_epsilon_n'],F(9,16));near(D(o['strict_tail_deltaSE'])**2,F(9,1024))
  passed.append(f'n{n} actual CLI matches all rational means and variance squares')
 rows=fixture(7)
 for eps,target,var in [(F(2,7),F(9,16),F(9,1024)),(F(3,7),F(3,8),F(9,256)),(F(2,7)-F(1,10**90),F(1),F(0)),(F(2,7)+F(1,10**90),F(9,16),F(9,1024)),(F(3,7)-F(1,10**90),F(9,16),F(9,1024)),(F(3,7)+F(1,10**90),F(3,8),F(9,256))]:
  o=run(rows,eps);near(o['strict_tail_d_gt_epsilon_n'],target);near(D(o['strict_tail_deltaSE'])**2,var)
 passed.append('n7 actual CLI rational threshold and both sides at 1e-90')
 for name,mut in [('extra_scaled_bin',lambda r:r[1]['tau2_distance_scaled_mean_L_hist'].append(0)),('nan_scaled_bin',lambda r:r[1]['tau2_distance_scaled_mean_L_hist'].__setitem__(2,float('nan'))),('upper_L2',lambda r:r[1]['tau2_distance_L2_groups'].__setitem__(0,100000))]:
  x=copy.deepcopy(rows);mut(x);run(x,'2/7',False);passed.append(name+' rejected')
receipt={'passed':passed,'aggregator_sha256':sha(AGG),'source_sha256':sha(SRC),'test_sha256':sha(pathlib.Path(__file__))}
ROOT.joinpath('hand-fixture-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
