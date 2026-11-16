#!/usr/bin/env python3
"""Strict Decimal-80 aggregation for the r022 tau=2 L-lift kernel."""
import argparse, hashlib, json, math, pathlib
from decimal import Decimal, getcontext
getcontext().prec=80
BASE=0xA022000000000000
# Rounded exp underflows only below half the least Double subnormal.
MINLOG=math.log(float.fromhex('0x0.0000000000001p-1022'))-math.log(2)
HALF_MIN_SUBNORMAL=Decimal(2) ** Decimal(-1075)
# Shortest-decimal parsing of Z plus one Double multiplication: four ulps of
# relative slack, with the required half-subnormal absolute rounding bound.
COMPLETE_REL_TOL=Decimal(2) ** Decimal(-48)
def dec(x): return Decimal(str(x))
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def bad(s): raise RuntimeError(s)
def finite(x):
    try:return dec(x).is_finite()
    except:return False
def ess(weights):
    s=sum(weights)
    return {'ESS':float(s*s/sum(x*x for x in weights)),'maxshare':float(max(weights)/s)}
def main():
 p=argparse.ArgumentParser();p.add_argument('raw');p.add_argument('--sidecar',required=True);p.add_argument('--source',required=True);a=p.parse_args()
 raw=pathlib.Path(a.raw);side=json.load(open(a.sidecar)); gotsha=sha(raw)
 if side.get('raw_sha256')!=gotsha or side.get('source_sha256')!=sha(a.source):bad('sidecar SHA mismatch')
 lines=[json.loads(x) for x in raw.read_text().splitlines() if x.strip()]
 if len(lines)<3 or lines[0].get('kind')!='metadata' or lines[-1].get('kind')!='run_end':bad('metadata/run_end ordering')
 m=lines[0]; rr=lines[1:-1]
 req=('implementation','base_source_sha256','n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps','lift_scope','lift_integer_bound')
 if any(k not in m for k in req) or m['implementation']!='r022-tau2-terminal-lift-complete-K2':bad('metadata schema')
 n,P,R=m['n'],m['particles'],m['replicates'];N=n*n
 expected={'n':n,'particles':P,'replicates':R,'tau_lo':2,'tau_hi':2,'domain':m['domain'],'method':m['method'],'move_steps':m['move_steps']}
 if side.get('expected_metadata')!=expected:bad('sidecar expected metadata mismatch')
 if not(n>=2 and n<=32 and P>0 and R>1 and m['tau_lo']==m['tau_hi']==2 and len(rr)==R):bad('metadata range/full-reps')
 if not finite(lines[-1].get('wall_seconds')):bad('run_end')
 if sorted(x.get('replicate',-1) for x in rr)!=list(range(R)):bad('replicate set')
 G=(P+255)//256; domain=int(m['domain'],16); seen=set(); vals=[]
 gf=('tau2_lift_L_groups','tau2_lift_L_distance2_groups','tau2_lift_L_F0_groups','tau2_lift_L_G0_groups','tau2_lift_group_LF_minus_2LG_minus_4LD')
 sf=('tau2_lift_scaled_mean_L','tau2_lift_scaled_mean_L_distance2','tau2_lift_scaled_mean_L_F0','tau2_lift_scaled_mean_L_G0')
 cf=('tau2_lift_complete_K2_Z','tau2_lift_complete_K2_D','tau2_lift_complete_K2_F','tau2_lift_complete_K2_G')
 core=('logZ','U_over_Z','D','F','G','F_over_Z','G_over_Z','F_minus_2G_minus_4D')
 for x in rr:
  if x.get('kind')!='replicate' or x.get('tau')!=2:bad('non-tau2 record')
  r=x['replicate']; want=(BASE+domain+n*1_000_000+2*10_000+r)&((1<<64)-1)
  if x.get('seed_hex','').lower().zfill(16)!=f'{want:016x}' or x['seed_hex'] in seen:bad('seed guard')
  seen.add(x['seed_hex'])
  if x.get('dead'):
   # On a free n-by-n grid with tau=2, the first root has a non-neighbour
   # available for the second root and every non-full prefix has boundary.
   # Hence a dead tau=2 trajectory is an implementation invariant failure,
   # not a valid zero contribution for this target.
   bad('impossible dead tau2 grid trajectory')
  if any(x.get(k) is None or not finite(x[k]) for k in core+sf+cf):bad('nondead numeric fields')
  if any(not isinstance(x.get(k),list) or len(x[k])!=G for k in gf):bad('group schema')
  if any(type(v) is not int for k in gf for v in x[k]):bad('non-integer group payload')
  zlog=float(x['logZ']); zdec=dec(x['logZ']).exp()
  if x['Z']==0 and not zlog<MINLOG:bad('false underflow')
  if x['Z']!=0 and abs(dec(x['Z'])-zdec)>dec('1e-12')*max(abs(dec(x['Z'])),abs(zdec))+HALF_MIN_SUBNORMAL:bad('Z/logZ')
  if abs(dec(x['F'])-2*dec(x['G'])-4*dec(x['D'])-dec(x['F_minus_2G_minus_4D']))>dec('1e-8')*max(dec(1),abs(dec(x['F'])),abs(2*dec(x['G'])),abs(4*dec(x['D']))):bad('base root identity')
  scaledres=dec(x['F_over_Z'])-2*dec(x['G_over_Z'])-4*dec(N)*dec(x['U_over_Z'])
  if abs(scaledres)>dec('1e-10')*max(dec(1),abs(dec(x['F_over_Z'])),abs(2*dec(x['G_over_Z'])),abs(4*dec(N)*dec(x['U_over_Z']))):bad('scaled base root identity')
  W=DD=FF=GG=0
  for g in range(G):
   cap=min(256,P-256*g);w,dd,ff,gg=(int(x[k][g]) for k in gf[:4]);res=int(x[gf[4]][g])
   if not(cap<=w<=cap*(N-3) and 0<=dd<=w*1922 and 0<=ff<=w*3844 and abs(gg)<=w*1922):bad('group bound')
   if ff-2*gg!=4*dd or res!=0:bad('group identity')
   W+=w;DD+=dd;FF+=ff;GG+=gg
  means=[dec(W)/P,dec(DD)/P,dec(FF)/P,dec(GG)/P]
  if any(abs(dec(x[k])-v)>dec('1e-12')*max(dec(1),abs(v)) for k,v in zip(sf,means)):bad('scaled means')
  complete=[zdec*v for v in means] # Decimal values used for actual aggregation.
  reported=[dec(x['Z'])*dec(v) for v in (x[k] for k in sf)]
  for k,v in zip(cf,reported):
   if x['Z']==0:
    if x[k]!=0:bad('underflowed complete not zero')
   elif abs(dec(x[k])-v)>COMPLETE_REL_TOL*abs(v)+HALF_MIN_SUBNORMAL:bad('complete lift field')
  vals.append((zdec,)+tuple(complete))
 if len(seen)!=R:bad('seed multiplicity')
 S=sum(v[1] for v in vals)
 if not S:bad('zero lift normalizer')
 Tz=sum(v[0] for v in vals);TD=sum(v[2] for v in vals);TF=sum(v[3] for v in vals);TG=sum(v[4] for v in vals)
 Q=TD/(N*S);V=TF/(8*S);C=TG/(4*S)
 def dse(component,theta,fac):
  h=[v[component]-fac*theta*v[1] for v in vals]
  # equivalent to z_r(d_r-fac*denom*l_r), without dividing an underflowed Z.
  return (dec(R)/dec(R-1)*sum((t/S)**2 for t in h)/(fac*fac)).sqrt()
 out={'kind':'tau2_lift_strict_aggregate','n':n,'P':P,'R':R,'tau2_normalizer_mean':str(Tz/R),'lift_full_normalizer_mean':str(S/R),'lift_full_D_mean':str(TD/R),'Q':str(Q),'V':str(V),'C':str(C),'deltaSE_Q':str(dse(2,Q,N)),'deltaSE_V':str(dse(3,V,8)),'deltaSE_C':str(dse(4,C,4)),'tau2_weight':ess([v[0] for v in vals]),'liftZ_weight':ess([v[1] for v in vals]),'note':'Finite-particle/SMC error and random-stream assumptions are not certified; ESS is diagnostic only.'}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
