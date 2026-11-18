"""Reconstruct underresolved boundary layers and test the extended cutoff."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from scipy.special import gammaln
import trace_scan as old
import trace_refined as new


def main():
    files=[json.loads(x.read_text()) for x in Path('trace').glob('shard-*.json')]
    index=sorted(((r['max_nested'],x['shard'],i) for x in files for i,r in enumerate(x['results'])),reverse=True)
    rows=[]
    for score,shard,rep in index[:10]:
        u,v,w,q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(750363+shard+15*rep)).random_base2(11).T
        p=2-np.exp2(-11+10.9*(shard+u)/15)
        c=np.exp2(-2.5+4.5*v)
        a=-1+4*w
        hlog=-(8+155*q)*np.log(2.)
        select=np.concatenate([np.arange(lo,min(2048,lo+384),41) for lo in range(0,2048,384)])
        m=[z[select] for z in [p,c,a,hlog]]
        y=old.normalized(*m,16);r=old.normalized(*m,32)
        i=select[np.argmax(abs(y-r))]
        assert abs(max(abs(y-r))-score)<1e-13
        model=[z[i:i+1] for z in [p,c,a,hlog]]
        vals=[float(new.normalized(*model,k)[0]) for k in (32,64,96)]
        rows.append(dict(profile=[float(x[0]) for x in model],old_max_disagreement=score,
                        refined=vals,max_fine_disagreement=abs(vals[-1]-vals[-2])))
    p=np.array([1.07,1.15,1.3,1.45,1.6,1.73])
    q0=1/p-.5
    c=np.full(len(p),1.2);a=np.full(len(p),.4)
    hlog=np.full(len(p),-112.)
    expected=np.exp(-q0*np.log(2*c)+gammaln(1-q0)+(2*q0)*hlog)
    old_ratio=old.normalized(p,c,a,hlog,32)/expected
    new_ratio=new.normalized(p,c,a,hlog,64)/expected
    out=dict(worst=rows,near_endpoint_fixed_cusp=dict(p=p.tolist(),hlog=-112,
                 original_ratio=old_ratio.tolist(),refined_ratio=new_ratio.tolist()))
    Path('refine_high.json').write_text(json.dumps(out,indent=2)+"\n")
    print(out['near_endpoint_fixed_cusp'])
    print(max(r['max_fine_disagreement'] for r in rows))

if __name__=='__main__':main()
