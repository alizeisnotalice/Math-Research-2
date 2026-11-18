"""Reconstruct highest nested-rule-discrepancy points and refine them locally."""
from pathlib import Path
import json

import numpy as np
from scipy.stats import qmc
from finite_packet import integral


def main():
    reports=[json.loads(p.read_text()) for p in Path("finite").glob("shard-*.json")]
    records=[(item["max_nested_quadrature_error"],report["shard"],k)
             for report in reports for k,item in enumerate(report["replicas"])]
    rows=[]
    for error,shard,repeat in sorted(records,reverse=True)[:5]:
        q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(795510+shard+15*repeat)).random_base2(16)
        ids=np.concatenate([start+np.arange(0,min(768,len(q)-start),40) for start in range(0,len(q),768)])
        u,v,_,_=q[ids].T
        d=np.exp(np.log(2.**-12)+(shard+u)/15*np.log(.75/(2.**-12)))
        lam=np.exp(np.log(16.)+v*np.log(64.))
        vals=[integral(lam,d,n) for n in (32,64,96,128)]
        arg=int(np.argmax(np.abs(vals[0]-vals[1])))
        assert abs(abs(vals[0][arg]-vals[1][arg])-error)<1e-11
        rows.append(dict(shard=shard,repeat=repeat,index=int(ids[arg]),delta=float(d[arg]),
                         lambda_=float(lam[arg]),order32_to64=error,
                         order64_to96=float(abs(vals[1][arg]-vals[2][arg])),
                         order96_to128=float(abs(vals[2][arg]-vals[3][arg])),
                         values=[str(v[arg]) for v in vals]))
    Path('worst-check.json').write_text(json.dumps(rows,indent=2)+"\n")
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    main()
