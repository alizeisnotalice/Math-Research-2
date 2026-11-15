#!/usr/bin/env python3
"""Small exact check of R46 final L1 half-deficit constants."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,os,sys

def main():
    assert os.path.ismount('/Volumes/data')
    ip=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json')
    data=json.loads(ip.read_text());m=sum(F(p[2]) for p in data['pts'])
    rows=[]
    for N in [6,7,10]:
        T=2*m*N*N;n=408*N*N;K=(1<<n)-1;loss=9*N-1
        upper_root=isqrt(T.numerator//T.denominator)+1
        assert F(upper_root)**2>T
        eta=F(loss,4*K*(upper_root+1))
        assert 0<eta<=1
        error=K*(4*eta*upper_root+4*eta*eta)
        assert error<=4*K*eta*(upper_root+1)==loss
        atomic_deficit=F(loss,m*N*N)
        assert atomic_deficit-error/T>=atomic_deficit/2
        rows.append({'N':N,'atom_count_upper':n,'mass':str(T),'sqrt_mass_rational_upper':upper_root,
                     'eta_formula':'(9*N-1)/(4*(2**n-1)*(sqrt_mass_rational_upper+1))',
                     'eta_denominator_bits':eta.denominator.bit_length(),
                     'area_error_at_most_9N_minus_1':True,'at_least_half_ratio_deficit_retained':True})
    out={'status':'PASS','python':sys.version,'single_process':True,'arithmetic':'Fraction/int/isqrt',
         'scope':'Three explicit rational choices support the general algebraic proof; no halo or smoothing file generated.',
         'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest(),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':rows}
    path=Path(__file__).parent/'results/reflection-edge-state-review/l1-half-gap.json'
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':len(rows),'max_atom_count_bound':max(r['atom_count_upper'] for r in rows),'output':str(path)}))

if __name__=='__main__':main()
