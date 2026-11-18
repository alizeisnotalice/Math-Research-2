"""Audit the protected-peak production residues and all merged upper counts."""
import hashlib
import json
from math import comb, factorial, prod
from pathlib import Path

from count_protected import independent_sizes, verified_primes


def sha(data):
    return hashlib.sha256(data).hexdigest()


def merged_integers(a, b, lists):
    f, t, v = lists[0]
    g, s, z = lists[1]
    N = a + b
    base = comb(N,a) * sum(f[i]*g[j] for i in range(3) for j in range(3-i))
    m1 = m2 = 0
    for i, j in ((0,2), (1,1), (2,0)):
        m1 += comb(N-1,a-1)*t[i]*g[j] + comb(N-1,a)*f[i]*s[j]
        m2 += comb(N-2,a-1)*((f[i]-t[i])*s[j] + t[i]*(g[j]-s[j]))
        m2 += comb(N-2,a-2)*v[i]*g[j] + comb(N-2,a)*f[i]*z[j]
    return base, m1, m2, base-m1-m2


def main():
    p = Path(__file__).resolve().parent
    s = json.loads((p/'residue-ledger.json').read_text())
    prime_path = p/'pratt-61bit.json'
    primes, nodes = verified_primes(prime_path)
    assert nodes == 40
    modulus = prod(primes)
    coefficients = [modulus//k * pow(modulus//k, -1, k) for k in primes]
    counts = {}
    jobs = {'local': 0, 'macbook': 0}
    root_source_hash = sha((p/'protected_mod.cpp').read_bytes())
    for name, batch in s.items():
        rec = batch['record']
        assert rec['primes'] == list(map(str,primes))
        assert rec['source_sha256'] == root_source_hash
        assert rec['prime_cert_sha256'] == sha(prime_path.read_bytes())
        assert rec['driver_sha256'] == sha((p/'count_protected.py').read_bytes())
        assert rec['regional_jobs'] == rec['completed_jobs'] == len(rec['receipts']) == len(batch['raw'])
        jobs[batch['machine']] += rec['completed_jobs']
        by_region = {}
        for row in rec['receipts']:
            h, ix = row['profile'], row['prime_index']
            raw = batch['raw'][f"{h}-p{ix}"].encode()
            assert sha(raw) == row['raw_sha256'] and row['exit_code'] == 0
            got = json.loads(raw)
            assert h == got['profile'] and ix == got['prime_index'] and int(got['prime']) == primes[ix]
            sizes = independent_sizes(''.join(str(int(w)-1) for w in h))
            assert got['states_by_layer'] == got['processed_by_layer'] == sizes
            assert got['independent_sets'] == sum(sizes)
            assert got['max_table_bytes'] <= 512*(1 << 20)
            by_region.setdefault(h,{})[ix] = got
        recovered = {}
        fields = ('peak_frequency_residues','unprotected_max_frequency_residues','bad_second_frequency_residues')
        for region, records in by_region.items():
            arrays = [[sum(c*x for c,x in zip(coefficients,col,strict=True)) % modulus
                       for col in zip(*(records[i][f] for i in range(4)),strict=True)]
                      for f in fields]
            f,t,v = arrays
            assert sum(f) == factorial(50)
            assert sum(t) == len(region)*factorial(49)
            assert all(0 <= v_i <= f_i-t_i <= f_i for f_i,t_i,v_i in zip(f,t,v,strict=True))
            assert arrays == [[int(x) for x in row] for row in rec['region_histograms'][region]]
            recovered[region] = arrays
        cut, = rec['upper_counts']
        d = name, cut['reflected_other_side']
        assert d[1] == ''.join(str(10-int(x)) for x in name)
        merged = merged_integers(50,50,[recovered[h] for h in d])
        keys = ('safe_peak_at_most_two_count','global_max_excluded_count',
                'global_second_excluded_count','A2_upper_count')
        assert merged == tuple(int(cut[k]) for k in keys)
        counts[name] = merged[-1]
        assert 0 <= counts[name] <= factorial(100)
    assert jobs == {'local':16,'macbook':20}
    assert min(counts,key=counts.get) == '5555555555'
    straight = json.loads((p/'straight-upper.json').read_text())
    assert int(straight['upper_counts'][0]['A2_upper_count']) == counts['5555555555']
    assert straight['source_sha256'] == root_source_hash
    print(json.dumps(dict(verified_jobs=jobs, primes=primes, pratt_nodes=nodes,
                          upper_counts=counts),indent=2))


if __name__ == '__main__':
    main()
