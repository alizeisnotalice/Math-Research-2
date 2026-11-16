"""Independent complete-log audit and finite stratified-ratio diagnostics."""
import json,hashlib,sys
from pathlib import Path
from decimal import Decimal,localcontext
p=Path(sys.argv[1]);raw=p.read_bytes();rows=[json.loads(s) for s in raw.splitlines()]
assert rows[0]['kind']=='metadata' and rows[-1]['kind']=='run_end'
assert sum(r['kind']=='metadata' for r in rows)==sum(r['kind']=='run_end' for r in rows)==1
m=rows[0];n=m['n'];N=n*n;R=m['replicates'];assert (n,m['particles'],R,m['tau_lo'],m['tau_hi'],m['method'])==(16,1048576,16,2,256,'strat')
rr=rows[1:-1];assert len(rr)==255*R
seen=set();states={};seedset=set();dead=[]
with localcontext() as ctx:
 ctx.prec=80
 D=lambda x:Decimal(str(x))
 for r in rr:
  assert r['kind']=='replicate';t,k=r['tau'],r['replicate'];key=t,k
  assert key not in seen and 2<=t<=N and 0<=k<R;seen.add(key)
  seed=(0xA022000000000000+int(m['domain'],16)+n*1000000+t*10000+k)&((1<<64)-1)
  assert int(r['seed_hex'],16)==seed and seed not in seedset;seedset.add(seed)
  assert r['scaled_U']==r['U_over_Z']
  if r['dead']:
   assert r['logZ'] is None and r['Z']==r['U']==r['U_over_Z']==0
   dead.append(key);z=u=D(0)
  else:
   ell=D(r['logZ']);q=D(r['U_over_Z']);assert ell.is_finite() and q.is_finite() and 0<=q<=D(2*(n-1)**2)/N
   z=ell.exp();u=z*q
   assert abs(D(r['Z'])-z)<=D('2e-12')*z
   assert abs(D(r['U'])-D(r['Z'])*q)<=D('2e-12')*z
  states[key]=(z,u)
 assert seen=={(t,k) for t in range(2,N+1) for k in range(R)}
 assert set(dead)=={(t,k) for t in [255,256] for k in range(R)}
 sz=sum(z for z,u in states.values());su=sum(u for z,u in states.values());Q=su/sz;Z=sz/R
 var=lambda xs:sum((x-sum(xs)/len(xs))**2 for x in xs)/(len(xs)-1)
 hs={t:[states[t,k][1]-Q*states[t,k][0] for k in range(R)] for t in range(2,N+1)}
 fixedse=(sum(var(v) for v in hs.values())/R).sqrt()/Z
 hz=[sum(hs[t][k] for t in hs) for k in range(R)];pairedse=(var(hz)/R).sqrt()/Z
 zrep=[sum(states[t,k][0] for t in hs) for k in range(R)]
 ess=sz**2/sum(z*z for z in zrep)
 wt={t:sum(states[t,k][0] for k in range(R))/sz for t in hs}
 out={'raw_sha256':hashlib.sha256(raw).hexdigest(),'complete_records':len(rr),'all_unique_seeds':len(seedset),'dead_records':len(dead),'Z_mean':str(Z),'Q_ratio_of_means':str(Q),'Q_delta_SE_fixed_strata':str(fixedse),'Q_delta_SE_paired_by_repeat':str(pairedse),'repeat_weight_ESS':str(ess),'repeat_max_weight':str(max(zrep)/sz),'estimated_second_birth_rank_mean':str(sum(t*w for t,w in wt.items())),'estimated_second_birth_tail':{str(k):str(sum(w for t,w in wt.items() if t>k)) for k in [2,3,4,5,8,16,32,64,128]},'note':'Finite descriptive diagnostics only. Fixed-strata delta SE uses independent-run sampling idealization; neither SE is a confidence bound or convergence/mixing certificate.'}
 print(json.dumps(out,indent=2))
