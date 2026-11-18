"""Exact trace replay of both host reductions and finite run diagnostics."""
from fractions import Fraction
from hashlib import sha256
from math import exp,isclose
from pathlib import Path
import argparse
import json

D=Path(__file__).resolve().parent


def read_batch(mode, traces=False):
    d=D/mode
    raw=d/'run.jsonl'
    r=json.loads((d/'finished.json').read_text())
    assert r['returncode']==r['caffeinate_returncode']==0 and r['waited_reaped']
    assert sha256(raw.read_bytes()).hexdigest()==r['raw_sha256']
    assert sha256((d/'resource.txt').read_bytes()).hexdigest()==r['resource_sha256']
    source=D/('contact.swift' if mode.startswith('contact') else 'transport.swift')
    if not source.exists():
        raise RuntimeError('generate contact.swift with make_contact.py before verifying')
    assert sha256(source.read_bytes()).hexdigest()==r['source_sha256']
    rows=[json.loads(s) for s in raw.read_text().splitlines()]
    m=rows[0]
    expected=r['parameters']
    assert (m['n'],m['particles'],m['replicates'],m['move_steps'],m['domain'])==tuple(
          expected[k] for k in ('n','P','R','L','domain'))
    assert len(rows)==2*m['replicates']+2==r['rows'] and rows[-1]['kind']=='run_end'
    for i in range(m['replicates']):
        rep,rb=rows[2*i+1:2*i+3]
        assert rep['replicate']==rb['replicate']==i and not rep['dead']
        seed=(0xA022000000000000+int(m['domain'],16)+m['n']*1_000_000+20000+i)%(1<<64)
        assert int(rep['seed_hex'],16)==seed
        h=rb['transport_histograms']
        assert len(h)==4 and h[2]==rb['mean_m_hist'][::-1]
        assert all(isclose(sum(k),rb['mean_m'],rel_tol=1e-12) for k in h)
    if traces:
        g=json.loads((d/'run.jsonl.trace.json').read_text())
        return rows,g['rows'][-1]
    return rows


def rho(p,q,n):
    return abs(p%n-q%n)+abs(p//n-q//n)


def contact_word(q,prefix,n):
    p=prefix[0]
    assert q==prefix[1]
    rest=[v for v in prefix if v!=q]
    idx=next(i for i,v in enumerate(rest) if rho(v,q,n)==1)
    K=rest[:idx]
    U={v for v in range(n*n) if all(rho(v,u,n)>1 for u in K)}
    assert q in U
    reached={q}
    while extra:=(U-reached)&{v for v in U if any(rho(v,u,n)==1 for u in reached)}:
        reached|=extra
    return idx+1,rho(p,q,n),max(rho(p,u,n) for u in reached),rest[idx] in U


def trace_gate(mode):
    rows,trace=read_batch(mode,traces=True)
    n=rows[0]['n']
    P=rows[0]['particles']
    D=2*n-2
    hist=[[Fraction(0) for _ in range(D+1)] for _ in range(4)]
    ranks=[0]*(n*n); joint=[0]*((D+1)**2); nonmerge=0
    r=trace['first_merge_rb_payload_ab_r1_r2']
    for i,word in enumerate(trace['prefix']):
        a,b,p,q=r[4*i:4*i+4]
        assert (p,q)==tuple(word[:2])
        wt=Fraction((a+b)*(a+b-1),2*a*b)
        x=p%n+p//n
        y=q%n+q//n
        nearest=min(x,D-x)+min(y,D-y)
        fixed=D-abs(x-y)
        adapt=D-rho(p,q,n)
        for v,at in zip(hist,(nearest,fixed,adapt,fixed-adapt),strict=True):
            v[at]+=wt
        j,d,far,grow=contact_word(q,word,n)
        ranks[j]+=j-1
        joint[d*(D+1)+far]+=j-1
        nonmerge+=int(grow)*(j-1)
    observed=rows[2]['transport_histograms']
    assert all(isclose(float(a/P),b,rel_tol=1e-12) for x,y in zip(hist,observed,strict=True)
               for a,b in zip(x,y,strict=True))
    if mode.startswith('contact'):
        obs=rows[2]['contact_counts']
        assert ranks==obs['contact_rank_L'] and joint==obs['contact_joint_source_max_distance_L']
        assert nonmerge==obs['first_contact_nonmerge_L'][0]
    return {'n':n,'particles':P,'lift_L_total':sum(ranks),'lift_L_nonmerge':nonmerge}


def legacy_equal(rows,other):
    assert rows[0]==other[0]
    for index in range(1,len(rows)-1,2):
        a,b=rows[index:index+2]
        x,y=other[index:index+2]
        assert all(a[k]==x[k] for k in ('seed_hex','logZ','tau2_lift_scaled_mean_L'))
        assert all(b[k]==y[k] for k in ('mean_m','mean_m_d2','mean_m_F0','mean_m_G0','mean_m_hist'))


def diagnostics(mode,rows):
    meta=rows[0]
    n=meta['n']
    base=max(x['logZ'] for x in rows[1:-1:2])
    weights=[exp(x['logZ']-base) for x in rows[1:-1:2]]
    rb=rows[2:-1:2]
    den=sum(u*x['mean_m'] for u,x in zip(weights,rb,strict=True))
    names=('old_nearest','fixed_distinct_corners','adaptive_distinct_corners','fixed_minus_adaptive')
    hist=[[sum(u*x['transport_histograms'][k][d] for u,x in zip(weights,rb,strict=True))/den
           for d in range(2*n-1)] for k in range(4)]
    averages={key:sum(d*x for d,x in enumerate(h)) for key,h in zip(names,hist,strict=True)}
    assert isclose(averages[names[1]]-averages[names[2]],averages[names[3]],abs_tol=1e-11)
    result=dict(n=n,P=meta['particles'],R=meta['replicates'],L=meta['move_steps'],
                domain=meta['domain'],mean_path_length=averages,finite_histograms=dict(zip(names,hist)),
                wall_seconds=rows[-1]['wall_seconds'],gpu_seconds=sum(x['gpu_seconds'] for x in rows[1:-1:2]))
    if mode.startswith('contact'):
        counts=[x['contact_counts'] for x in rb]
        size=2*n-1
        matrices=[x['contact_joint_source_max_distance_L'] for x in counts]
        normal=sum(w*sum(z) for w,z in zip(weights,matrices,strict=True))
        target_lines=[]
        for t in (n//2+1,n-1,2*n-2):
            num=sum(w*sum(z[d*size+f] for d in range(t) for f in range(t,size)) for w,z in zip(weights,matrices))
            eligible=sum(w*sum(z[:t*size]) for w,z in zip(weights,matrices))
            assert 0<=num<=eligible<=normal
            target_lines.append({'target_at_least':t,'near_source_fraction':eligible/normal,
                                 'reachable_in_near_sources':num/eligible})
        result['rank_lift_contact']={'targets':target_lines,
               'nonmerge_fraction':sum(w*x['first_contact_nonmerge_L'][0] for w,x in zip(weights,counts))/normal}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--gate',choices=('gate','contact-gate'))
    p.add_argument('--r052',type=Path)
    args=p.parse_args()
    if args.gate:
        print(json.dumps(trace_gate(args.gate),indent=2))
    else:
        assert args.r052 is not None
        prev=[json.loads(x) for x in args.r052.read_text().splitlines()]
        current={mode:read_batch(mode) for mode in ('replay','n24','contact-n16')}
        legacy_equal(current['replay'],prev)
        legacy_equal(current['contact-n16'],prev)
        print(json.dumps({'gates':{k:trace_gate(k) for k in ('gate','contact-gate')},
                          'identical_earlier_legacy_replays':2,
                          'diagnostics':{k:diagnostics(k,v) for k,v in current.items()}},indent=2))
