"""Exact odd-size four-quadrant numerator and rectangular root gate."""
from fractions import Fraction
from functools import cache
from math import factorial, prod
import hashlib
import json
from pathlib import Path
import platform

P = Path(__file__).resolve().parent

def roots(height, width):
    N = height*width
    adj = [sum(1 << v for v in range(N)
               if abs(v%width-u%width)+abs(v//width-u//width)==1) for u in range(N)]

    @cache
    def count(S):
        if S == 0:
            return factorial(N)
        cover = S
        total = 0
        for b in bits(S):
            cover |= adj[b.bit_length()-1]
            total += count(S^b)
        answer, remainder = divmod(total,cover.bit_count())
        assert remainder == 0
        return answer

    answer = [0]*N

    def visit(S,candidates):
        if S:
            value = count(S) if S.bit_count()%2 else -count(S)
            for b in bits(S):
                answer[b.bit_length()-1] += value
        while candidates:
            b = candidates & -candidates
            candidates ^= b
            visit(S|b,candidates&~adj[b.bit_length()-1])

    visit(0,(1<<N)-1)
    return answer,count.cache_info().currsize

def bits(mask):
    while mask:
        b=mask&-mask
        yield b
        mask^=b

def sha(content):
    return hashlib.sha256(content).hexdigest()

def main():
    source = P.parents[1]/'r046/layered'
    receipt = json.loads((P/'modular-batch.json').read_text())
    old = receipt['summary']
    assert sha((source/'root_mod.cpp').read_bytes()) == old['counter_source_sha256']
    assert sha((source/'pratt-61bit.json').read_bytes()) == old['prime_cert_sha256']
    full, independent_sets = roots(6,5)
    assert independent_sets == sum(old['planned']['555555']['states_by_layer'])
    primes = list(map(int,old['primes']))
    modulus = prod(primes)
    assert modulus > factorial(30)
    crt = [modulus//p * pow(modulus//p,-1,p) for p in primes]
    rows = []
    for rec in receipt['receipts']:
        i = rec['prime_index']
        raw = receipt['stdout'][f'555555-p{i}'].encode()
        assert sha(raw)==rec['stdout_sha256']
        assert rec['exit_code']==0
        values=json.loads(raw)
        assert values['processed_by_layer']==old['planned']['555555']['states_by_layer']
        assert int(values['prime'])==primes[i]
        rows.append((i,list(map(int,values['root_residues']))))
    assert len(rows)==3
    residues = [row for _,row in sorted(rows)]
    independent_crt = [sum(a*b for a,b in zip(col,crt,strict=True))%modulus
                       for col in zip(*residues,strict=True)]
    assert independent_crt == list(map(int,old['root_counts']['555555']))
    checked=0
    for v,c in enumerate(independent_crt):
        if v%6<2:
            assert c==full[5*(v//6)+v%6] and c>0
            checked += 1
        else:
            assert c==0
    assert checked==12
    a5=json.loads((P.parents[1]/'r048/n10-gluing/n10-quadrants.json').read_text())['all_root_counts']
    a6=list(map(int,json.loads((P.parents[2]/'computation/results/independent-peaks-n6-r004.json').read_text())['single_peak_root_counts']))
    left=sum(a5[x]*full[x] for x in range(2))
    right=sum(full[5*x]*a6[x] for x in range(4))
    fact=prod(factorial(k) for k in (25,30,30,36))
    total=factorial(121)//fact*left*right
    answer=dict(python=platform.python_version(),n=11,distance_lower=6,independent_sets=independent_sets,
                rectangle_root_counts=full,verified_safe_roots=checked,root_products=[left,right],
                constructed_words=total,uniform_unconditional=str(Fraction(left*right,fact)))
    print(json.dumps(answer,indent=2))

if __name__=='__main__':
    main()
