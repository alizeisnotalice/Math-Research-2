#!/usr/bin/env python3
"""Log-scaled Decimal aggregation for completed r023 tau=2 distance runs.

All ratios are evaluated after replacing each terminal normalizer Z_r by
exp(logZ_r-max_s logZ_s).  This preserves every ratio and delta-SE while
avoiding a dependence on IEEE-double normalizer underflow at large n.
"""
import argparse, hashlib, json, math, pathlib
from fractions import Fraction
from decimal import Decimal, getcontext
# More guard digits than the historical Decimal-80 reducer are used because
# translating all log weights by a large common constant changes only the
# rounding path, never the intended ratio.
getcontext().prec=110
BASE=0xA022000000000000
# Rounded exp underflows only below half the least Double subnormal.
MINLOG=math.log(float.fromhex('0x0.0000000000001p-1022'))-math.log(2)
HALF_MIN_SUBNORMAL=Decimal(2) ** Decimal(-1075)
# Shortest-decimal parsing of Z plus one Double multiplication: four ulps of
# relative slack, with the required half-subnormal absolute rounding bound.
COMPLETE_REL_TOL=Decimal(2) ** Decimal(-48)
def dec(x): return Decimal(str(x))
def epsfrac(s): return Fraction(s)  # exact for a/b and finite decimal text
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def bad(s): raise RuntimeError(s)
def finite(x):
    try:return dec(x).is_finite()
    except:return False
def ess(weights):
    s=sum(weights)
    return {'ESS':float(s*s/sum(x*x for x in weights)),'maxshare':float(max(weights)/s)}
def main():
 p=argparse.ArgumentParser();p.add_argument('raw');p.add_argument('--sidecar',required=True);p.add_argument('--source',required=True);p.add_argument('--epsilon',type=str);a=p.parse_args()
 raw=pathlib.Path(a.raw);side=json.load(open(a.sidecar)); gotsha=sha(raw)
 if side.get('raw_sha256')!=gotsha or side.get('source_sha256')!=sha(a.source):bad('sidecar SHA mismatch')
 lines=[json.loads(x) for x in raw.read_text().splitlines() if x.strip()]
 if len(lines)<3 or lines[0].get('kind')!='metadata' or lines[-1].get('kind')!='run_end':bad('metadata/run_end ordering')
 m=lines[0]; rr=lines[1:-1]
 req=('implementation','base_source_sha256','n','particles','replicates','tau_lo','tau_hi','domain','method','move_steps','lift_scope','lift_integer_bound')
 if any(k not in m for k in req) or m['implementation']!='r023-tau2-terminal-L-distance-histogram' or 'distance_scope' not in m:bad('metadata schema')
 n,P,R=m['n'],m['particles'],m['replicates'];N=n*n
 expected={'n':n,'particles':P,'replicates':R,'tau_lo':2,'tau_hi':2,'domain':m['domain'],'method':m['method'],'move_steps':m['move_steps']}
 if side.get('expected_metadata')!=expected:bad('sidecar expected metadata mismatch')
 if not(n>=2 and n<=32 and P>0 and R>1 and m['tau_lo']==m['tau_hi']==2 and len(rr)==R):bad('metadata range/full-reps')
 if not finite(lines[-1].get('wall_seconds')):bad('run_end')
 if sorted(x.get('replicate',-1) for x in rr)!=list(range(R)):bad('replicate set')
 G=(P+255)//256; domain=int(m['domain'],16); seen=set(); vals=[]; dhist=[]; dl2=[]; Mdist=2*n-1
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
  # The producer accumulates logZ and multiplicative Z along different
  # Double paths.  Actual n=16 records differ by about 6.1e-15 relatively,
  # so this is a deliberately conservative *consistency* guard, not an
  # assertion that logZ is a correctly rounded logarithm of Z.
  if x['Z']!=0 and abs(dec(x['Z'])-zdec)>dec('1e-12')*max(abs(dec(x['Z'])),abs(zdec))+HALF_MIN_SUBNORMAL:bad('Z/logZ')
  # Raw root moments are separately emitted in the producer.  Unlike the
  # aggregation weights, they are formed with its reported Double Z.  Check
  # that complete tuple as well, so a self-consistent but unrelated D/F/G
  # payload cannot pass merely through F-2G=4D.
  if x['Z']==0:
   if any(dec(x[k])!=0 for k in ('D','F','G','F_minus_2G_minus_4D')):bad('underflowed root tuple not zero')
  else:
   rz=dec(x['Z']); targets=(rz*dec(N)*dec(x['U_over_Z']),rz*dec(x['F_over_Z']),rz*dec(x['G_over_Z']))
   for k,v in zip(('D','F','G'),targets):
    if abs(dec(x[k])-v)>COMPLETE_REL_TOL*max(abs(v),abs(dec(x[k])))+HALF_MIN_SUBNORMAL:bad('raw root scale tuple')
  if abs(dec(x['F'])-2*dec(x['G'])-4*dec(x['D'])-dec(x['F_minus_2G_minus_4D']))>dec('1e-8')*max(dec(1),abs(dec(x['F'])),abs(2*dec(x['G'])),abs(4*dec(x['D']))):bad('base root identity')
  scaledres=dec(x['F_over_Z'])-2*dec(x['G_over_Z'])-4*dec(N)*dec(x['U_over_Z'])
  if abs(scaledres)>dec('1e-10')*max(dec(1),abs(dec(x['F_over_Z'])),abs(2*dec(x['G_over_Z'])),abs(4*dec(N)*dec(x['U_over_Z']))):bad('scaled base root identity')
  W=DD=FF=GG=0
  for g in range(G):
   cap=min(256,P-256*g);w,dd,ff,gg=(int(x[k][g]) for k in gf[:4]);res=int(x[gf[4]][g])
   if not(cap<=w<=cap*(N-3) and 0<=dd<=w*1922 and 0<=ff<=w*3844 and abs(gg)<=w*1922):bad('group bound')
   if ff-2*gg!=4*dd or res!=0:bad('group identity')
   W+=w;DD+=dd;FF+=ff;GG+=gg
  hg=x.get('tau2_distance_L_hist_groups');l2g=x.get('tau2_distance_L2_groups')
  if not isinstance(hg,list) or not isinstance(l2g,list) or len(hg)!=G or len(l2g)!=G:bad('distance group schema')
  H=[];L2=0
  for g in range(G):
   cap=min(256,P-256*g); row=hg[g]
   if not isinstance(row,list) or len(row)!=Mdist or any(type(v) is not int for v in row) or type(l2g[g]) is not int:bad('distance integer payload')
   if row[0] or row[1] or sum(row)!=int(x['tau2_lift_L_groups'][g]):bad('distance hist identity')
   if any(v<0 or v>cap*(N-3) for v in row):bad('distance hist bound')
   q2=int(l2g[g]);w=int(x['tau2_lift_L_groups'][g])
   if not (w<=q2<=(N-3)*w and q2<=cap*(N-3)*(N-3)):bad('L2 bound')
   # integer Cauchy: sum L_i^2 >= ceil(W^2/cap)
   if q2 < (w*w+cap-1)//cap:bad('L2 Cauchy')
   H.append(row);L2+=q2
  flat=[sum(H[g][d] for g in range(G)) for d in range(Mdist)]
  sh=x.get('tau2_distance_scaled_mean_L_hist')
  if not isinstance(sh,list) or len(sh)!=Mdist or any(not finite(v) for v in sh) or not finite(x.get('tau2_distance_scaled_mean_L2')):bad('distance scaled schema')
  if any(abs(dec(sh[d])-dec(flat[d])/P)>dec('1e-12') for d in range(Mdist)) or abs(dec(x['tau2_distance_scaled_mean_L2'])-dec(L2)/P)>dec('1e-12'):bad('distance scaled fields')
  dhist.append(flat);dl2.append(L2)
  means=[dec(W)/P,dec(DD)/P,dec(FF)/P,dec(GG)/P]
  if any(abs(dec(x[k])-v)>dec('1e-12')*max(dec(1),abs(v)) for k,v in zip(sf,means)):bad('scaled means')
  complete=[zdec*v for v in means] # Only for checking reported Double fields.
  reported=[dec(x['Z'])*dec(v) for v in (x[k] for k in sf)]
  for k,v in zip(cf,reported):
   if x['Z']==0:
    if x[k]!=0:bad('underflowed complete not zero')
   elif abs(dec(x[k])-v)>COMPLETE_REL_TOL*abs(v)+HALF_MIN_SUBNORMAL:bad('complete lift field')
  # Delay exponentiation until every record has been audited.  The retained
  # payload is independent of the absolute normalizer scale.
  vals.append((dec(x['logZ']),)+tuple(means))
 if len(seen)!=R:bad('seed multiplicity')
 # Decimal has a broad exponent range, but this explicit common shift is the
 # object actually used below.  It makes normalizer underflow impossible in
 # the ratio calculation and permits a reproducible log-normalizer report.
 maxlog=max(v[0] for v in vals)
 vals=[( (v[0]-maxlog).exp(),)+tuple((v[0]-maxlog).exp()*x for x in v[1:]) for v in vals]
 S=sum(v[1] for v in vals)
 if not S:bad('zero lift normalizer')
 Tz=sum(v[0] for v in vals);TD=sum(v[2] for v in vals);TF=sum(v[3] for v in vals);TG=sum(v[4] for v in vals)
 Q=TD/(N*S);V=TF/(8*S);C=TG/(4*S)
 def dse(component,theta,fac):
  h=[v[component]-fac*theta*v[1] for v in vals]
  # equivalent to z_r(d_r-fac*denom*l_r), without dividing an underflowed Z.
  return (dec(R)/dec(R-1)*sum((t/S)**2 for t in h)/(fac*fac)).sqrt()
 histfull=[sum(vals[r][0]*dec(dhist[r][d])/P for r in range(R)) for d in range(Mdist)]
 probs=[h/S for h in histfull]
 hse=[]
 for d,theta in enumerate(probs): hse.append((dec(R)/dec(R-1)*sum(((vals[r][0]*dec(dhist[r][d])/P-theta*vals[r][1])/S)**2 for r in range(R))).sqrt())
 E2=(sum(vals[r][0]*(dec(dl2[r])/P+3*vals[r][1]/vals[r][0]) for r in range(R)))/(2*S)
 ese=(dec(R)/dec(R-1)*sum(((vals[r][0]*(dec(dl2[r])/P+3*vals[r][1]/vals[r][0])-2*E2*vals[r][1])/(2*S))**2 for r in range(R))).sqrt()
 tail=None;tailse=None
 if a.epsilon is not None:
  ep=epsfrac(a.epsilon); inds=[d for d in range(Mdist) if d*ep.denominator>n*ep.numerator]
  tail=sum(probs[d] for d in inds); tr=[vals[r][0]*sum(dec(dhist[r][d]) for d in inds)/P for r in range(R)]
  tailse=(dec(R)/dec(R-1)*sum(((tr[r]-tail*vals[r][1])/S)**2 for r in range(R))).sqrt()
 out={'kind':'tau2_distance_logscaled_aggregate','n':n,'P':P,'R':R,
      'normalizer_log_shift':str(maxlog),
      'log_tau2_normalizer_mean':str(maxlog+(Tz/R).ln()),
      'log_lift_full_normalizer_mean':str(maxlog+(S/R).ln()),
      'tau2_normalizer_mean_scaled_by_exp_minus_shift':str(Tz/R),
      'lift_full_normalizer_mean_scaled_by_exp_minus_shift':str(S/R),
      'Q':str(Q),'V':str(V),'C':str(C),'deltaSE_Q':str(dse(2,Q,N)),
      'distance_L_weighted_probabilities':[str(x) for x in probs],'distance_deltaSE':[str(x) for x in hse],
      'E2_tau':str(E2),'deltaSE_E2_tau':str(ese),'epsilon':a.epsilon,
      'strict_tail_d_gt_epsilon_n':None if tail is None else str(tail),
      'strict_tail_deltaSE':None if tailse is None else str(tailse),
      'tau2_weight':ess([v[0] for v in vals]),'liftZ_weight':ess([v[1] for v in vals]),
      'note':'Finite diagnostic only; all ratios use Decimal exp(logZ-maxlog); strict tail uses exact Fraction cross multiplication.'}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
