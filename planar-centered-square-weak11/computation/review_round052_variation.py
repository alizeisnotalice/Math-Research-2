#!/usr/bin/env python3
"""Author replay and exact scalar checks for the R52 geometry/variation bridge.

Does not repeat 45,633,320 bbox scan. Pins producer code, replays full endpoint
union certificate, then independently bounds its ten-term derivative and tests
rational contained radii. Default paths are source evidence in isolated clones;
command-line args allow archived paths to be supplied for reproduction.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import importlib.util
import json
from math import isqrt
import os
from pathlib import Path
import platform
import sys
import time
sys.dont_write_bytecode=True
PRODUCER_SHA='37bbe43363c673f9118b34fbeaf2301710acf3ff70e3fb3d2046b74cad863960'
BASE=Path('/Volumes/data/CodexProjects/T2/.local/workspaces')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root(q,bits=160):
    D=1<<bits;n=isqrt((q.numerator*D*D)//q.denominator)
    low=F(n,D);high=low if low*low==q else F(n+1,D)
    assert low*low<=q<=high*high
    return low,high

def scalar_derivative(terms,vlo,vhi):
    lo=hi=F(0)
    for row in terms:
        q,b,c=map(F,(row['q0'],row['q1'],row['coefficient']))
        if not b or not c:continue
        ql,qh=sorted((q+b*vlo,q+b*vhi));assert ql>0
        rl,rh=root(ql)[0],root(qh)[1];k=c*b/2
        a,z=k/rh,k/rl;lo+=min(a,z);hi+=max(a,z)
    return lo,hi

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--producer',type=Path,default=BASE/'pcw-kernel-audit/planar-centered-square-weak11/computation/position_mass_endpoint_union.py')
    ap.add_argument('--slab',type=Path,default=BASE/'pcw-position-slab-impl/planar-centered-square-weak11/computation/round052/position-mass-slab.compact.json.gz')
    ap.add_argument('--certificate',type=Path,default=BASE/'pcw-kernel-audit/planar-centered-square-weak11/computation/results/round052-endpoint-union/certificate.json')
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();start=time.process_time()
    assert sha(args.producer)==PRODUCER_SHA
    slab=json.loads(gzip.decompress(args.slab.read_bytes()));cert=json.loads(args.certificate.read_text())
    spec=importlib.util.spec_from_file_location('r052_reviewed_producer',args.producer);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    replay=mod.serial(mod.certify(slab,160))
    checkkeys=[k for k in replay if k not in ('cpu_seconds',)]
    for k in checkkeys:assert replay[k]==cert[k],k
    T=F(cert['dt_radius']);center=scalar_derivative(cert['beta_root_terms'],F(0),F(0));cell=scalar_derivative(cert['beta_root_terms'],-T,T)
    assert center==tuple(map(F,cert['beta_t_center_bounds']))
    assert cell==tuple(map(F,cert['beta_t_cell_bounds']))
    assert cell[1]<0 and F(cert['mass_t_derivative'])==0
    M,lam=F(cert['mass']),F(cert['lambda'])
    exact_mixed=tuple(lam*x/M for x in cell)
    assert exact_mixed==tuple(map(F,cert['mixed_J_z_t_cell_bounds']))
    assert F(1,50_000_000)<F(cert['dz_radius'])
    assert F(3,5000)<F(cert['dt_radius'])
    # Outward rounded decimal fractions are deliberately coarser than the
    # exact interval, and these inequalities certify their displayed use.
    display_lower=F('-6.121211246656e-9');display_upper=F('-6.121211242977e-9')
    assert display_lower<exact_mixed[0]<=exact_mixed[1]<display_upper<0
    # A local maximum of a(t)+z*b(t) needs b'=0: for any nonzero b',
    # det Hessian = -b'^2 < 0, independently of a''+z*b''.
    hessian_checks=0
    for a in [F(-10),F(-1),F(0),F(7,2)]:
        for b in [F(-3),F(-1,17),F(1,19),F(4)]:
            assert a*0-b*b<0;hessian_checks+=1
    # The necessary condition must not be overstated as b identically zero:
    # -t^2+z*t^2 <= 0 for |z|<=1/2 and all tested rational t.
    nonflat_max_checks=0
    for z in [F(-1,2),F(0),F(1,2)]:
        for t in [F(-1,10),F(0),F(1,10)]:
            assert -t*t+z*t*t<=0;nonflat_max_checks+=1
    out={'status':'passed','scope':'author replay of complete joint endpoint/symbolic-union stage; independent exact root-sum scalar derivative and normalization checks; no fresh bbox enumeration',
       'producer_sha256':sha(args.producer),'slab_sha256':sha(args.slab),'certificate_sha256':sha(args.certificate),
       'replayed_fields':checkkeys,'replayed_fields_count':len(checkkeys),'box_count':cert['box_count'],
       'independent_beta_t_center_bounds':list(map(str,center)),'independent_beta_t_cell_bounds':list(map(str,cell)),
       'independent_mixed_cell_bounds':list(map(str,exact_mixed)),
       'contained_exact_dz_radius':'1/50000000','contained_exact_dt_radius':'3/5000',
       'outward_decimal_mixed_bounds':[str(display_lower),str(display_upper)],
       'hessian_determinant_examples':hessian_checks,'nonflat_local_max_examples':nonflat_max_checks,
       'python':sys.version,'platform':platform.platform(),'cpu_seconds':time.process_time()-start,
       'script_sha256':sha(Path(__file__)),'actual_command':[sys.executable,*sys.argv]}
    assert os.path.ismount('/Volumes/data')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','replayed_fields_count','box_count','cpu_seconds','script_sha256')},indent=2))
if __name__=='__main__':main()
