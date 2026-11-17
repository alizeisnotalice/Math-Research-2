#!/usr/bin/env python3
"""Independent integer CRT and recursive primality-certificate verifier."""
import argparse,hashlib,json,math
from pathlib import Path

ap=argparse.ArgumentParser();ap.add_argument("certificate",type=Path);ap.add_argument("residues",type=Path,nargs=4);args=ap.parse_args()
cert=json.loads(args.certificate.read_text());nodes=cert["nodes"];done=set()
def prime(p):
    if p in done:return
    node=nodes[str(p)];assert int(node["p"])==p
    if p==2:
        assert node=={"p":"2","base":"two"};done.add(p);return
    factors={int(q):int(e) for q,e in node["p_minus_1_factorization"].items()}
    assert factors and all(2<=q<p and e>=1 for q,e in factors.items())
    assert math.prod(q**e for q,e in factors.items())==p-1
    for q in factors:prime(q)
    a=int(node["witness"]);assert 1<a<p and pow(a,p-1,p)==1
    assert all(math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in factors)
    done.add(p)
primes=list(map(int,cert["roots"]));assert len(primes)==len(set(primes))==4
for p in primes:prime(p)
# Independently enumerate the cardinality distribution by row transfer.
valid=[r for r in range(256) if not r&(r<<1)]
dp={0:[1]+[0]*24}
for _ in range(6):
    new={r:[0]*25 for r in valid}
    for r in valid:
        shift=r.bit_count()
        for prev,poly in dp.items():
            if not r&prev:
                for k in range(25-shift):new[r][k+shift]+=poly[k]
    dp=new
hist=[sum(poly[k] for poly in dp.values()) for k in range(25)]
count,peak=sum(hist),max(hist)
assert count==851302029 and peak==162959262
rows=[json.loads(path.read_text()) for path in args.residues]
rows=sorted(rows,key=lambda x:x["prime_index"])
assert [r["prime_index"] for r in rows]==list(range(4))
for i,r in enumerate(rows):
    assert r["n"]==8 and r["algorithm"]=="modular-layer-rolling-reenumerated-hash"
    assert int(r["prime"])==primes[i]>56 and 0<=int(r["residue"])<primes[i]
    assert r["independent_sets"]==count and r["max_layer_count"]==peak
    assert r["max_modeled_payload_bytes"]<=r["payload_limit_bytes"]==7*2**30
res=[int(r["residue"]) for r in rows]
x,M=0,1
for r,p in zip(res,primes):
    x+=M*((r-x)*pow(M,-1,p)%p);M*=p
alternate=sum(r*(M//p)*pow(M//p,-1,p) for r,p in zip(res,primes))%M
assert x==alternate and all(x%p==r for r,p in zip(res,primes))
upper=math.factorial(56);assert M>upper and 0<x<=upper
out={"verified":True,"scope":"Exact finite boundary completion count for n=8; no full A2(8) count or asymptotic claim", "n":8,"H_toprow":str(x),"CRT_modulus":str(M),"upper_bound_56_factorial":str(upper),"CRT_unique_within_proved_bound":True,"primality_certificate_nodes_verified":len(done),"independent_set_layer_histogram":hist,"residues":[str(r) for r in res],"primes":[str(p) for p in primes],"input_sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in [args.certificate,*args.residues]}}
print(json.dumps(out,indent=2))
