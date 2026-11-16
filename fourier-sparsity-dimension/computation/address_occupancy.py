#!/usr/bin/env python3
"""Exact spectral-coset occupancy/energy checks for address functions."""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
from exact_baseline import basis


def run(k):
    t=1 << k
    # Frequency mask: k selector bits followed by t independent data bits.
    spectrum={a | (1 << (k+i)): (-1)**((a&i).bit_count()%2)
              for i in range(t) for a in range(t)}
    assert len(spectrum)==t*t and len(basis(spectrum))==t+k
    records=[]
    for j in range(k+1):
        K=1 << j
        wmask=K-1
        cosets=defaultdict(list)
        for a,c in spectrum.items():
            cosets[a & ~wmask].append((a & wmask,c))
        ell=len(cosets)
        assert ell==t*t//K
        active=[0]*K
        energies=[0]*K
        for terms in cosets.values():
            assert len(terms)==K
            vals=[sum(c*(-1)**((a&b).bit_count()%2) for a,c in terms) for b in range(K)]
            p=Fraction(sum(v!=0 for v in vals),K)
            mu=Fraction(sum(v*v for v in vals),K*t*t)
            fourth=Fraction(sum(v**4 for v in vals),K*t**4)
            assert p==Fraction(1,K) and mu==Fraction(K,t*t)
            assert p==max(mu,Fraction(1,len(terms)))
            assert p==mu*mu/fourth
            for b,v in enumerate(vals):
                active[b]+=int(v!=0)
                energies[b]+=v*v
        assert all(v==t*t//(K*K) for v in active)
        assert all(v==t*t for v in energies)
        assert Fraction(sum(active),K)==Fraction(ell*ell,t*t)
        records.append({'j':j,'query_rank':j,'cosets':ell,'fiber_size':K,
                        'active_per_branch':active[0],'p':str(Fraction(1,K)),
                        'mu':str(Fraction(K,t*t)),
                        'energy_uncertainty_equality':True,'second_fourth_moment_equality':True})
    return {'selector_bits':k,'s':t*t,'d':t+k,'stages':records}


ap=argparse.ArgumentParser();ap.add_argument('--max-k',type=int,default=6)
ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
r={'schema':'a28-address-occupancy-v1','python':sys.version,'platform':platform.platform(),
   'arithmetic':'integer numerators, denominator t; Fraction for moments',
   'seed':None,'tolerance':0,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   'cases':[run(k) for k in range(a.max_k+1)]}
a.out.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'status':'passed','k_range':[0,a.max_k],'stages':sum(len(c['stages']) for c in r['cases'])}))
