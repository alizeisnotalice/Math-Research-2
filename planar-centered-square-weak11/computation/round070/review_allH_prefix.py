#!/usr/bin/env python3
"""Independent exhaustive subset vs tied-prefix audit with exact arithmetic."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
import itertools
import json
import os
from pathlib import Path
import platform
import random
import time


def digest(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def subset_le(r,M,m):
    # r <= sqrt(M+m)-sqrt(m)/2 iff q>=0 and q²>=r²m.
    q=M+3*m/4-r*r
    return q>=0 and q*q>=r*r*m


def sqrt_interval(v,bits):
    # Rational enclosure proven by integer isqrt; no floating arithmetic.
    den=1<<bits
    k=isqrt((v.numerator*den*den)//v.denominator)
    lo=F(k,den)
    return lo, lo if lo*lo==v else F(k+1,den)


def prefix_le(r,M,m,stats):
    # Independent adaptive comparison using enclosures of the original radicals.
    # Exact equality recognized first; doubling otherwise must separate the roots.
    q=M+3*m/4-r*r
    equal=q>=0 and q*q==r*r*m
    if equal:
        stats['radical_equalities']+=1
        return True
    bits=4
    while True:
        slo,shi=sqrt_interval(M+m,bits)
        tlo,thi=sqrt_interval(m,bits)
        lower=slo-thi/2
        upper=shi-tlo/2
        stats['max_isolation_bits']=max(stats['max_isolation_bits'],bits)
        if r<=lower:
            return True
        if r>upper:
            return False
        bits*=2


def audit_case(points,masses,U,m,x,stats):
    L=sum(masses)
    remaining=[i for i in range(len(points)) if not U>>i&1]
    distances={i:max(abs(x[j]-points[i][j]) for j in (0,1)) for i in remaining}
    exhaustive=[]
    for mask in range(1,1<<len(remaining)):
        H=[remaining[j] for j in range(len(remaining)) if mask>>j&1]
        M=sum(masses[i] for i in H)
        r=max(distances[i] for i in H)
        bad=M+m>L and subset_le(r,M,m)
        exhaustive.append((H,M,r,bad))
    prefixes=[]
    for r in sorted(set(distances.values())):
        H=[i for i in remaining if distances[i]<=r]
        M=sum(masses[i] for i in H)
        geom=prefix_le(r,M,m,stats)
        assert geom==subset_le(r,M,m)
        prefixes.append((H,M,r,M+m>L and geom))
    assert any(row[3] for row in exhaustive)==any(row[3] for row in prefixes)
    # Stronger check: every individual violating H expands to its tied prefix.
    by_r={row[2]:row for row in prefixes}
    for H,M,r,bad in exhaustive:
        if bad:
            P,MP,rP,badP=by_r[r]
            assert set(H)<=set(P) and MP>=M and rP==r and badP
    stats['subset_rows']+=len(exhaustive)
    stats['prefix_rows']+=len(prefixes)
    stats['cases_with_violation']+=any(row[3] for row in exhaustive)
    stats['tied_cases']+=len(prefixes)<len(remaining)
    return not any(row[3] for row in exhaustive)


def main():
    assert __debug__
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--balanced-input',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args();start=time.perf_counter()
    data=json.loads(Path(args.input).read_text())
    stats={'subset_rows':0,'prefix_rows':0,'cases_with_violation':0,'tied_cases':0,
           'radical_equalities':0,'max_isolation_bits':0}
    scalar=[]
    for row in data['fixed_scalar_cases']:
        r,M,m=[F(row[k]) for k in ('r','M','m')]
        a=subset_le(r,M,m);b=prefix_le(r,M,m,stats)
        q=M+3*m/4-r*r
        assert a==b==row['comparison_le']
        assert (q>=0 and q*q==r*r*m)==row['equality']
        scalar.append({**row,'PASS':True})
    triggers=[]
    for row in data['fixed_trigger_cases']:
        M,m,L,r=[F(row[k]) for k in ('M','m','L','r')]
        assert subset_le(r,M,m)
        assert (M+m>L)==row['trigger']
        triggers.append({**row,'PASS':True})
    rng=random.Random(data['seed']);random_rows=[]
    for case in range(data['random_cases']):
        n=rng.randint(*data['rank_range'])
        raw=rng.sample(list(itertools.product(range(-8,9),repeat=2)),n)
        points=[tuple(F(v,data['coordinate_denominator']) for v in p) for p in raw]
        masses=[F(rng.randint(1,16),data['mass_denominator']) for _ in range(n)]
        U=rng.randint(1,(1<<n)-2)
        a=sum(masses[i] for i in range(n) if U>>i&1);L=sum(masses)
        m=a+(L-a)*F(rng.randint(1,9),10)
        x=tuple(F(rng.randint(-12,12),8) for _ in (0,1))
        before=(stats['subset_rows'],stats['prefix_rows'])
        passed=audit_case(points,masses,U,m,x,stats)
        random_rows.append({'case':case,'points':[[str(v) for v in p] for p in points],
            'masses':list(map(str,masses)),'U_mask':U,'m':str(m),'x':list(map(str,x)),
            'subsets_checked':stats['subset_rows']-before[0],
            'prefixes_checked':stats['prefix_rows']-before[1],'compatible':passed})
    # An explicit full-tie case. The radius1/2 prefix includes both tied points.
    assert not audit_case([(F(0),F(0)),(F(1,2),F(0)),(-F(1,2),F(0)),(F(0),F(1,4))],
               [F(1,4)]*4,1,F(1,2),(F(0),F(0)),stats)
    # Test actual necessity, separately from random algebraic carrier fixtures.
    bdata=json.loads(Path(args.balanced_input).read_text())
    pp=[tuple(map(F,p)) for p in bdata['points']];ww=list(map(F,bdata['initial_masses']))
    Smask=bdata['S_mask'];L=sum(ww[i] for i in range(len(pp)) if Smask>>i&1)
    eligible=[]
    for mask in range(1,1<<len(pp)):
        ids=[i for i in range(len(pp)) if mask>>i&1]
        mass=sum(ww[i] for i in ids)
        span=max(max(pp[i][j] for i in ids)-min(pp[i][j] for i in ids) for j in (0,1))
        if mass>span*span:
            eligible.append((mask,ids,mass,span))
    maxima=[mask for mask,ids,mass,span in eligible if mass==max(r[2] for r in eligible)]
    assert maxima==[Smask] and L==1
    Sids=[i for i in range(len(pp)) if Smask>>i&1]
    spoints=[pp[i] for i in Sids];smasses=[ww[i] for i in Sids]
    real_labels=[];receivers=0
    for mask,ids,mass,span in eligible:
        if not mask&Smask or not mask&~Smask:
            continue
        U=sum(1<<j for j,i in enumerate(Sids) if mask>>i&1)
        assert 0<U<(1<<len(Sids))-1 and mass<=L
        # Find dyadic r strictly between span/2 and sqrt(mass)/2.
        bits=2
        while True:
            lo,hi=sqrt_interval(mass,bits)
            r=(span+lo)/4
            if r>span/2 and 4*r*r<mass:
                break
            bits*=2
        rect=[(max(pp[i][j] for i in ids)-r,min(pp[i][j] for i in ids)+r) for j in (0,1)]
        points_to_check=list(itertools.product(*[( (3*a+b)/4,(a+b)/2,(a+3*b)/4) for a,b in rect]))
        for x in points_to_check:
            assert all(4*max(abs(x[j]-pp[i][j]) for j in (0,1))**2<mass for i in ids)
            assert audit_case(spoints,smasses,U,mass,x,stats)
            receivers+=1
        real_labels.append(mask)
    assert len(eligible)==13 and len(real_labels)==4 and receivers==36
    output={'status':'PASS','scope':'independent all-subset/tied-prefix equivalence and actual maximum-block necessity checks; not W area',
        'python':platform.python_version(),'pid':os.getpid(),'source_sha256':digest(__file__),
        'input_sha256':digest(args.input),'balanced_input_sha256':digest(args.balanced_input),
        'seed':data['seed'],'random_cases':data['random_cases'],'stats':stats,
        'scalar_cases':scalar,'trigger_cases':triggers,'random_fixtures':random_rows,
        'actual_balanced':{'all_subsets':(1<<len(pp))-1,'eligible_count':len(eligible),
            'maximum_masks':maxima,'mixed_labels_checked':real_labels,'strict_rational_receivers':receivers},
        'arithmetic':'Fraction throughout; prefix radical order by adaptive exact dyadic enclosures with exact equality handling; no fixed precision or float acceptance',
        'elapsed_seconds':time.perf_counter()-start}
    Path(args.output).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','stats':stats,'actual_receivers':receivers,
                      'elapsed_seconds':output['elapsed_seconds']}))

if __name__=='__main__': main()
