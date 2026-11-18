"""Independent CRT and output-integrity check for the first n=9 cut."""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

cert=Path(__file__).resolve().parent
sys.path.insert(0,str(cert/'layered'))
from run_crt import independent_sizes, safe_roots, verified_primes

data=cert/'n9'
summary=json.loads((data/'summary.json').read_text())
read_sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
launch=json.loads((data/'launch.json').read_text())
for name,expected in launch['sha256'].items():
    assert read_sha(cert/'layered'/name)==expected
assert read_sha(data/'summary.json')==json.loads((data/'transfer-receipt.json').read_text())['summary_sha256']
primes,_=verified_primes(cert/'layered/pratt-61bit.json')
P=math.prod(primes)
assert P>math.factorial(41)
pres=[(P//p)*pow(P//p,-1,p) for p in primes]
vectors={}
rss=[]
keys=set()
for line in (data/'receipts.jsonl').read_text().splitlines():
    rec=json.loads(line)
    h,i=rec['profile'],rec['prime_index']
    assert (h,i) not in keys and rec['exit_code']==0
    keys.add((h,i))
    path=data/f'{h}-p{i}.json'
    assert read_sha(path)==rec['stdout_sha256']
    assert read_sha(path.with_suffix('.resource'))==rec['resource_sha256']
    value=json.loads(path.read_text())
    assert value['processed_by_layer']==value['states_by_layer']==independent_sizes(h)
    assert int(value['prime'])==primes[i] and value['profile']==h
    vectors[h,i]=list(map(int,value['root_residues']))
    rss.append(rec['max_rss_bytes'])
totals={}
for h in ('654343456','345656543'):
    assert len(vectors[h,0])==81
    area=sum(map(int,h))
    recovered=[sum(u*v for u,v in zip(values,pres))%P for values in zip(*(vectors[h,i] for i in range(3)))]
    assert recovered==list(map(int,summary['root_counts'][h]))
    assert {i for i,c in enumerate(recovered) if c}==safe_roots(h)
    assert sum(recovered)<=math.factorial(area)
    totals[h]=sum(recovered)
assert len(keys)==6
cut=math.comb(81,40)*math.prod(totals.values())
assert all(int(value)==cut for value in summary['individual_cut_counts'].values())
assert summary['max_active_table_payload_bytes']==3072*2**20
print(json.dumps(dict(count=cut,unconditional_probability_exact=str(Fraction(cut,math.factorial(81))),
                     regional_counts=totals,workers_finished=len(keys),max_rss=max(rss)),indent=2))
