"""Complete-stratum root moments; Decimal arithmetic, descriptive errors only."""
from collections import defaultdict
from decimal import Decimal as D, getcontext
from pathlib import Path
import hashlib
import json
import sys

getcontext().prec = 80

def near(x, y, tol=D('5e-12')):
    return x == y if not x or not y else abs(x-y) <= tol * max(abs(x), abs(y))

def summarize(path):
    path = Path(path)
    data = [json.loads(line, parse_float=D) for line in path.read_text().splitlines()]
    meta, end = data[0], data[-1]
    assert meta['kind'] == 'metadata' and end['kind'] == 'run_end', 'Incomplete run'
    assert meta['implementation'] == 'r020-soa-gpu-root-moments'
    n = meta['n']; N = n*n; R = meta['replicates']; domain = int(meta['domain'], 16)
    assert meta['tau_lo'] == 2 and meta['tau_hi'] == N and R > 1
    rows = data[1:-1]
    assert len(rows) == (N-1)*R and all(row['kind'] == 'replicate' for row in rows)
    groups = defaultdict(list); keys = set(); dead = 0; underflow = 0
    for row in rows:
        t, rep = row['tau'], row['replicate']
        assert 2 <= t <= N and 0 <= rep < R and (t, rep) not in keys
        keys.add((t, rep))
        seed = (0xA022000000000000 + domain + n*1000000 + t*10000 + rep) % (1 << 64)
        assert int(row['seed_hex'], 16) == seed
        if row['dead']:
            dead += 1
            assert row['Z'] == 0 and row['logZ'] is None
            assert all(row[k] is None for k in ['U','D','F','G','U_over_Z','F_over_Z','G_over_Z','V','C'])
            groups[t].append({k: D(0) for k in ['Z','U','D','F','G']})
            continue
        assert row['logZ'] is not None and row['Z'] >= 0
        z = D(row['Z']); q = D(row['U_over_Z']); fz = D(row['F_over_Z']); gz = D(row['G_over_Z'])
        assert D(0) <= q <= D(2)*(n-1)**2/N
        assert near(fz-2*gz, 4*N*q) and near(D(row['V']), fz/8) and near(D(row['C']), gz/4)
        if z == 0:
            underflow += 1
            z = D(row['logZ']).exp()
        else:
            assert near(z, D(row['logZ']).exp(), D('5e-10'))
            for key, expected in [('U',z*q),('D',z*N*q),('F',z*fz),('G',z*gz)]:
                assert near(D(row[key]), expected), (t, rep, key)
        groups[t].append({'Z':z, 'U':z*q, 'D':z*N*q, 'F':z*fz, 'G':z*gz})
    assert set(groups) == set(range(2,N+1)) and all(len(g)==R for g in groups.values())
    sums = {key: sum((sum(x[key] for x in g)/R for g in groups.values()), D(0)) for key in ['Z','U','D','F','G']}
    assert sums['Z'] > 0 and near(sums['F']-2*sums['G'],4*sums['D'])
    ratios = {'Q':sums['U']/sums['Z'], 'V_over_n2':sums['F']/(8*N*sums['Z']), 'C_over_n2':sums['G']/(4*N*sums['Z'])}
    errors = {}
    for name, key, scale in [('Q','U',1),('V_over_n2','F',8*N),('C_over_n2','G',4*N)]:
        ratio = ratios[name]*scale
        variance = D(0)
        for g in groups.values():
            h = [x[key]-ratio*x['Z'] for x in g]; avg = sum(h)/R
            variance += sum((x-avg)**2 for x in h)/(R*(R-1))
        errors[name] = str(variance.sqrt()/(scale*sums['Z']))
    weights = []
    for t,g in sorted(groups.items()):
        z = [x['Z'] for x in g]; total = sum(z)
        weights.append({'tau':t,'sum_Z':str(total),'ESS':str(total**2/sum(x*x for x in z)) if total else None,
                        'max_share':str(max(z)/total) if total else None})
    return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'metadata':meta,
            'coverage_records':len(rows),'dead':dead,'live_underflows_reconstructed_from_logZ':underflow,
            'ratios':{k:str(v) for k,v in ratios.items()},'sum_stratum_means':{k:str(v) for k,v in sums.items()},
            'descriptive_delta_SE':errors,'stratum_weight_diagnostics':weights,
            'scope':'Finite numerical FK ratio-of-means diagnostics. No certified confidence interval or mixing guarantee.'}

if __name__ == '__main__':
    print(json.dumps(summarize(sys.argv[1]), indent=2, default=str))
