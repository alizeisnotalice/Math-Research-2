"""Audit nearly-full-cap inverse conditioning by extended precision."""
import json
from pathlib import Path
import numpy as np
import mpmath as mp
from thresholds import get_cases,evaluate
import thresholds_refined as new


def main():
    mp.mp.dps=58
    records=[json.loads(p.read_text()) for p in Path('thresholds').glob('shard-*.json')]
    index=sorted(((r['max_absolute_threshold_error'],w['shard'],i)
                  for w in records for i,r in enumerate(w['results'])),reverse=True)
    out=[]
    for score,s,i in index[:8]:
        l,r,m,mu,a0,b0,p0=get_cases(s,i,11)
        a,b,p,_,est=evaluate(l,r,m,mu)
        j=int(np.argmax(np.maximum(abs(a-a0),abs(b-b0))))
        assert abs(max(abs(a[j]-a0[j]),abs(b[j]-b0[j]))-score)<1e-14
        ll,rr,mm,ff=map(lambda v:mp.mpf(float(v[j])),(l,r,m,mu))
        phi=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
        cdf=lambda x:mp.erfc(-x/mp.sqrt(2))/2
        c0=cdf(ll);c1=cdf(-rr)
        z=lambda p:mp.sqrt(2)*mp.erfinv(2*p-1)
        def thresholds(P):
            return z(c0+P),-z(c1+(mm-P))
        def mom(P):
            a,b=thresholds(P)
            return phi(ll)-phi(a)+phi(b)-phi(rr)
        assert mom(mm)<ff<mom(mp.mpf(0))
        P,Q=mp.mpf(0),mm
        for _ in range(140):
            mid=(P+Q)/2
            if mom(mid)>ff:P=mid
            else:Q=mid
        aa,bb=thresholds((P+Q)/2)
        vv=mm+ll*phi(ll)-aa*phi(aa)+bb*phi(bb)-rr*phi(rr)
        vals=new.evaluate(l[j:j+1],r[j:j+1],m[j:j+1],mu[j:j+1])
        out.append(dict(l=float(ll),r=float(rr),mass=float(mm),first=float(ff),
                        cap_share=float(mm/(cdf(rr)-cdf(ll))),base_absolute_threshold_disagreement=score,
                        revised_absolute_threshold_error_vs_precise_target=float(max(abs(vals[0][0]-aa),abs(vals[1][0]-bb))),
                        base_absolute_threshold_error_vs_precise_target=float(max(abs(a[j]-aa),abs(b[j]-bb))),
                        target_rounding_shift_vs_original_left=float(abs(aa-a0[j])),
                        base_second_absolute_error_vs_precise_target=float(abs(est[j]-vv)),
                        precise_thresholds_for_rounded_target=[str(aa),str(bb)]))
    Path('conditioning.json').write_text(json.dumps(out,indent=2)+"\n")
    print(out[0])
if __name__=='__main__':main()
