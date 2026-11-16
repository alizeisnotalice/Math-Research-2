#!/usr/bin/env python3
"""Decimal-80 log-sum-exp summary, counting dead repeats as exact zero."""
from decimal import Decimal, getcontext
import json, sys, hashlib
getcontext().prec=80
for path in sys.argv[1:]:
    rows=[json.loads(s, parse_float=Decimal) for s in open(path)]; reps=[x for x in rows if x.get('kind')=='replicate'];R=len(reps)
    assert rows[0]['kind']=='metadata' and rows[-1]['kind']=='run_end'
    assert rows[0]['tau_lo']==rows[0]['tau_hi'], 'this audit is for one stratum only'
    assert R==rows[0]['replicates'] and sorted(x['replicate'] for x in reps)==list(range(R))
    for x in reps:
        expected=(0xA022000000000000+int(rows[0]['domain'],16)+rows[0]['n']*1000000+x['tau']*10000+x['replicate'])&((1<<64)-1)
        assert x['tau']==rows[0]['tau_lo'] and int(x['seed_hex'],16)==expected
        assert (x['logZ'] is None)==x['dead']
    live=[x for x in reps if not x['dead']]
    if live:
        logs=[Decimal(str(x['logZ'])) for x in live];m=max(logs);w=[(x-m).exp() for x in logs]
        total=(m.exp())*sum(w);mean=total/Decimal(R);logmean=(m+(sum(w)/Decimal(R)).ln())
        z=[x.exp() for x in logs]+[Decimal(0)]*(R-len(live))
        shares=[x/sum(w) for x in w];ess=Decimal(1)/sum(x*x for x in shares);maxshare=max(shares)
        logmin,logmax=min(logs),max(logs)
        qs=[Decimal(str(x["U_over_Z"])) for x in live if x.get("U_over_Z") is not None]
        weighted_q=(sum(x.exp()*Decimal(str(y["U_over_Z"])) for x,y in zip(logs,live) if y.get("U_over_Z") is not None)/sum(x.exp() for x,y in zip(logs,live) if y.get("U_over_Z") is not None)) if qs else None
        unweighted_q=sum(qs)/Decimal(len(qs)) if qs else None
    else:
        mean=Decimal(0);logmean=None;z=[Decimal(0)]*R;ess=Decimal(0);maxshare=Decimal(0);logmin=logmax=None;weighted_q=unweighted_q=None
    variance=sum((x-mean)**2 for x in z)/Decimal(R-1) if R>1 else Decimal(0)
    se=(variance/Decimal(R)).sqrt() if R else Decimal(0)
    out={'file':path.rsplit('/',1)[-1],'n':rows[0]['n'],'P':rows[0]['particles'],'R':R,'tau_lo':rows[0]['tau_lo'],'tau_hi':rows[0]['tau_hi'],'domain':rows[0]['domain'],'dead':R-len(live),'log_mean_Z':None if logmean is None else str(logmean),'mean_Z_decimal':str(mean),'sample_variance_Z_decimal':str(variance),'descriptive_SE_Z_decimal':str(se),'repeat_logZ_ESS':str(ess),'repeat_max_share':str(maxshare),'logZ_min':None if logmin is None else str(logmin),'logZ_max':None if logmax is None else str(logmax),'unweighted_mean_U_over_Z':None if unweighted_q is None else str(unweighted_q),'Z_weighted_U_over_Z_single_tau_only':None if weighted_q is None else str(weighted_q),'wall_seconds':rows[-1]['wall_seconds'],'gpu_seconds_sum':sum(x['gpu_seconds'] for x in reps)}
    out['sha256']=hashlib.sha256(open(path,'rb').read()).hexdigest();out['scope']='single-stratum finite numerical diagnostic; no mixing or coverage certificate'
    print(json.dumps(out,sort_keys=True,default=str))
