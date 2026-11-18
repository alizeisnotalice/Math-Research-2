"""Crosscheck five-modulus regional histograms and the strict original-law bound."""
import hashlib
import json
from math import comb, factorial, prod
from pathlib import Path
from fractions import Fraction

from count_protected import independent_sizes, merge, verified_primes

P = Path(__file__).resolve().parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def independent_merge(a,b,left,right):
    f,t,v=left
    g,s,z=right
    n=a+b
    raw=comb(n,a)*sum(f[i]*g[j] for i in range(3) for j in range(3-i))
    max_exclude=0
    second_exclude=0
    for i,j in ((0,2),(1,1),(2,0)):
        max_exclude += comb(n-1,a-1)*t[i]*g[j]+comb(n-1,a)*f[i]*s[j]
        second_exclude += comb(n-2,a-2)*v[i]*g[j]+comb(n-2,a)*f[i]*z[j]
        second_exclude += comb(n-2,a-1)*((f[i]-t[i])*s[j]+t[i]*(g[j]-s[j]))
    return raw,max_exclude,second_exclude,raw-max_exclude-second_exclude

def main():
    s=json.loads((P/'complete-batch.json').read_text())
    primes,nodes=verified_primes(P/'pratt-61bit.json')
    assert nodes==50
    modulus=prod(primes)
    fractions=[modulus//p*pow(modulus//p,-1,p) for p in primes]
    h,k='55555566666','66666655555'
    assert modulus>factorial(61)
    assert len(s)==10
    by_region={h:{},k:{}}
    fields=('peak_frequency_residues','unprotected_max_frequency_residues','bad_second_frequency_residues')
    for item in s:
        r=item['receipt']
        meta=r['receipt']
        assert r['outcome']=='complete' and meta['exit_code']==0
        profile,index=meta['profile'],meta['prime_index']
        assert item['machine']==('local' if profile==h else 'macbook')
        assert sha(item['residue_json'].encode())==meta['raw_sha256']
        assert sha(item['resource'].encode())==meta['resource_sha256']
        assert sha((P/'protected_mod.cpp').read_bytes())==r['source_sha256']
        data=json.loads(item['residue_json'])
        assert data['processed_by_layer']==data['states_by_layer']==r['planned']['states_by_layer']
        expected=independent_sizes(''.join(str(int(c)-1) for c in profile))
        assert expected==data['states_by_layer']
        assert int(data['prime'])==primes[index] and data['profile']==profile
        assert index not in by_region[profile]
        by_region[profile][index]=data
    histograms={}
    for profile,residues in by_region.items():
        assert set(residues)==set(range(5))
        m=sum(map(int,profile))
        arrays=[[sum(x*y for x,y in zip(v,fractions,strict=True))%modulus
                 for v in zip(*(residues[i][field] for i in range(5)),strict=True)]
                for field in fields]
        f,t,v=arrays
        assert sum(f)==factorial(m) and sum(t)==11*factorial(m-1)
        assert all(0<=j<=a-b<=a for a,b,j in zip(f,t,v,strict=True))
        histograms[profile]=arrays
    a,b,c,d=independent_merge(60,61,histograms[h],histograms[k])
    data=merge(h,histograms)
    keys=('safe_peak_at_most_two_count','global_max_excluded_count','global_second_excluded_count','A2_upper_count')
    assert (a,b,c,d)==tuple(int(data[key]) for key in keys)
    numerator=json.loads((P.parent/'odd-quadrants/n11-numerator.json').read_text())['constructed_words']
    lower=Fraction(numerator,d)
    print(json.dumps(dict(n=11,pratt_nodes=nodes,five_moduli=list(primes),recovered_histograms=histograms,
               valid_workers=len(s),A2_upper=d,denominator_upper_over_factorial=str(Fraction(d,factorial(121))),
               uniform_far_threshold='11/2',far_conditional_lower=str(lower),
               A2_upper_decimal=d/factorial(121),lower_decimal=float(lower)),indent=2))

if __name__=='__main__':
    main()
