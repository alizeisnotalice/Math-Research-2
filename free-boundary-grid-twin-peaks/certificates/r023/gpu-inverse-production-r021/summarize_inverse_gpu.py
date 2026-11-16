#!/usr/bin/env python3
"""Full one-peak GPU batches, integer group checks, Decimal weighted ratios."""
import argparse,json,hashlib,math
from pathlib import Path
from decimal import Decimal as D,localcontext

def summarize(path):
    data=[json.loads(s,parse_float=D) for s in path.read_text().splitlines()]
    meta=data[0]; assert meta['kind']=='metadata'
    n,P,R,L=(int(meta[k]) for k in ('n','particles','replicates','move_steps'))
    N=n*n; assert meta['tau_lo']==meta['tau_hi']==N+1 and R>1
    assert data[-1]['kind']=='run_end'
    rows=[r for r in data if r['kind']=='replicate']
    assert len(rows)==R and {r['replicate'] for r in rows}==set(range(R))
    domain=int(meta['domain'],16)
    def close(x,y):
        assert abs(x-y)<=D('5e-12')*max(abs(x),abs(y)),(x,y)
    with localcontext() as c:
        c.prec=80
        zs=[]; ms={2:[],3:[]}; underflow=0
        for r in rows:
            assert r['tau']==N+1 and not r['dead']
            seed=(0xA022000000000000+domain+n*1000000+(N+1)*10000+r['replicate'])%(1<<64)
            assert int(r['seed_hex'],16)==seed
            z=D(r['Z']); assert z>=0 and z.is_finite()
            if not z:
                z=D(r['logZ']).exp();underflow+=1
            zs.append(z)
            assert r['K2'] is None
            for d in (2,3):
                groups=r[f'inverse_m{d}_groups'];assert len(groups)==(P+255)//256
                assert all(isinstance(g,int) and g>=0 for g in groups)
                bound=4*math.prod(range(N-d+1,N))
                for i,g in enumerate(groups): assert g<=min(256,P-i*256)*bound
                mean=D(sum(groups))/P
                close(mean,D(r[f'mean_m{d}']))
                close(mean,D(r[f'scaled_mean_m{d}']))
                if r['Z']:close(z*mean,D(r[f'Zmean_m{d}']))
                else:assert r[f'Zmean_m{d}']==0
                ms[d].append(z*mean)
        zsum=sum(zs); zm=zsum/R
        result={'input':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'n':n,'P':P,'R':R,'L':L,'domain_hex':meta['domain'],'mean_Z_onepeak_probability':str(zm),'mean_Z_times_factorial_A1':str(zm*math.factorial(N)),'nondead_underflows_restored':underflow,'replicate_weight_ESS':str(zsum*zsum/sum(z*z for z in zs)),'largest_replicate_weight_share':str(max(zs)/zsum),'ratios':{}}
        for d in (2,3):
            ratio=sum(ms[d])/zsum
            se=(sum((m-ratio*z)**2 for m,z in zip(ms[d],zs))/D(R*(R-1))).sqrt()/zm
            result['ratios'][str(d)]={'A2_distance_over_A1':str(ratio),'delta_SE':str(se)}
        result['scope']='Finite numerical SMC diagnostic; ratio bias and mixing are not certified.'
        return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=summarize(a.input);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
