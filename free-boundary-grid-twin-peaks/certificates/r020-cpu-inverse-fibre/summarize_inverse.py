#!/usr/bin/env python3
"""Validate complete one-peak SMC batches; Decimal ratio-of-means diagnostics.
SE is a replicate delta-method diagnostic, not a certified confidence interval.
"""
import argparse, csv, hashlib, json
from decimal import Decimal as D, localcontext
from pathlib import Path


def summarize(path):
    meta = json.loads(Path(str(path)+'.metadata.json').read_text())
    with path.open() as f:
        rows = list(csv.DictReader(f))
    n, P, R, L = (int(meta[k]) for k in ('n','P','R','L'))
    assert n >= 2 and P > 0 and R > 1 and L >= 0
    assert len(rows) == R and {int(r['rep']) for r in rows} == set(range(R))
    tau = n*n+1
    seed_base = int(meta['seed_base'])
    def close(x,y):
        assert abs(x-y) <= D('5e-14')*max(abs(x),abs(y)), (x,y)
    with localcontext() as ctx:
        ctx.prec=80
        for r in rows:
            assert int(r['n'])==n and int(r['tau'])==tau and int(r['L'])==L
            seed = (seed_base+tau*1000003+int(r['rep'])*1009+L*1000000007) % (1<<64)
            assert int(r['seed'])==seed
            assert 0 <= int(r['accept']) <= int(r['proposals'])
            # At rank N the terminal vertex is unique. No terminal MH.
            assert int(r['proposals']) == P*L*(n*n-1)
            assert 0 <= int(r['dead'])
            z=D(r['Z'])
            assert z.is_finite() and z>=0
            for d in (2,3):
                m=D(r[f'Zmean_m{d}'])
                assert m.is_finite() and m>=0
                if z:
                    close(m/z,D(r[f'weighted_ratio_m{d}']))
                else:
                    assert m==0
        zs=[D(r['Z']) for r in rows]
        total=sum(zs); assert total>0
        mean=total/R
        z_se=(sum((z-mean)**2 for z in zs)/D(R*(R-1))).sqrt()
        result={'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'metadata':meta, 'replicates':R, 'mean_Z':str(mean),'mean_Z_delta_SE':str(z_se),
            'replicate_weight_ESS':str(total**2/sum(z*z for z in zs)),
            'largest_replicate_weight_share':str(max(zs)/total),
            'zero_normalizer_replicates':sum(z==0 for z in zs),
            'dead_particle_rank_events':sum(int(r['dead']) for r in rows),
            'minimum_particle_parent_ESS':min(float(r['min_parent_ess']) for r in rows),
            'proposals':sum(int(r['proposals']) for r in rows),
            'accepted':sum(int(r['accept']) for r in rows), 'moments':{}}
        for d in (2,3):
            ms=[D(r[f'Zmean_m{d}']) for r in rows]
            ratio=sum(ms)/total
            residual=[m-ratio*z for m,z in zip(ms,zs)]
            # Sample mean residual is zero by construction.
            se=(sum(x*x for x in residual)/D(R*(R-1))).sqrt()/mean
            result['moments'][str(d)]={'ratio_of_means':str(ratio),'delta_SE':str(se)}
        result['scope']='Finite SMC diagnostic only; no certified confidence interval, mixing bound or asymptotic conclusion.'
        return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('csv',type=Path,nargs='+'); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); out={'runs':[summarize(path) for path in a.csv]}
    a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
