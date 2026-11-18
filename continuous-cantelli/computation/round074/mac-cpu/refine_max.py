"""Reconstruct the worst nested order disagreements and refine those cases."""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.stats import qmc
import mpmath as mp
import cantor_gibbs as cg
from critical_audit import renewal

def main():
    rec=[json.loads(p.read_text()) for p in Path("cantor").glob("shard-*.json")]
    index=sorted(((r['max_nested'],x['shard'],i) for x in rec for i,r in enumerate(x['results'])),reverse=True)
    for order in (128,192,256):
        n,w=roots_legendre(order)
        cg.RULES[order]=((n+1)/2,w/2)
    results=[]
    for score,shard,rep in index[:6]:
        u,v,w,q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(682413+shard+15*rep)).random_base2(8).T
        b=3+np.floor(10*u).astype(int)
        k=2+np.floor(23*v).astype(int)
        a=.45+1.55*(shard+w)/15
        a[::5]=2
        c=np.exp2(-3+6*q)
        sel=np.concatenate([np.arange(lo,min(256,lo+22),7) for lo in range(0,256,22)])
        parameters=[z[sel] for z in (b,k,a,c)]
        x=cg.coefficient(*parameters,48)[0];y=cg.coefficient(*parameters,80)[0]
        win=int(np.argmax(abs(x-y)))
        ix=sel[win]
        assert abs(abs(x[win]-y[win])-score)<3e-13
        model=[z[ix:ix+1] for z in (b,k,a,c)]
        vals=[float(cg.coefficient(*model,j)[0][0]) for j in (48,80,128,192,256)]
        row=dict(shard=shard,repeat=rep,index=int(ix),base=int(b[ix]),k=int(k[ix]),alpha=float(a[ix]),a=float(c[ix]),
                 rule_orders=[48,80,128,192,256],coefficients=vals,
                 high_order_gap=abs(vals[-1]-vals[-2]))
        if a[ix]==2:
            mpval=renewal(int(b[ix]),mp.mpf(float(c[ix])),k=int(k[ix]))[2]
            row['independent_exact_gaussian']=str(mpval)
            row['fine_absolute_error']=float(abs(vals[-1]-mpval))
        results.append(row)
    Path('refine_max.json').write_text(json.dumps(results,indent=2)+"\n")
    print(results)

if __name__=='__main__': main()
