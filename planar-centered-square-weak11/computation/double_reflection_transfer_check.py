#!/usr/bin/env python3
"""R54 independent exact tail and threshold-gap transfer check; no N scan."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,sys

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert os.path.ismount('/Volumes/data')
    base=Path(__file__).resolve().parents[1]
    source=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json')
    d=json.loads(source.read_text()); px,py=F(d['Px']),F(d['Py'])
    m=sum(F(p[2]) for p in d['pts']);rho=m/(px*py)
    short=rho*(1+px/15)*(1+py/30)
    long=rho*F(3,2)*(1+py/(6*px))
    c=max(short,long); eps=F(1,1024)
    smoothed=c*(1+eps/15)**2
    assert short<1 and long<1 and 6*px>30
    assert smoothed<1
    # Every N>=7 has T=2W>12Px. These fixed guard counts cover
    # both center halo expansion and all r<15 source dependencies.
    assert 2*(12*px)>15+eps and 4*py>15+eps
    cert=base/'computation/results/double-reflection-closure/exact-check.json'
    old=json.loads(cert.read_text())
    assert F(old['short_R15_tail_bound'])==short
    assert F(old['long_radius_tail_bound'])==long
    p=tuple(map(F,old['existing_R41_partial_coefficient']['bx_plus_true_gammaR_interval']))
    threshold=(-p[1],-p[0])
    assert threshold[1]<F('3.46882944')
    assert threshold[0]>F('3.46882943')
    out={
        'status':'PASS','scope':'Independent exact tails and sufficient smoothing/guard constants; no new seam area, N scan, or full weak-norm bound.',
        'python':sys.version,'arithmetic':'fractions.Fraction','random_seed':None,'floating_tolerance':None,
        'R':'15','N_domain':'all integers N>=7','Px':str(px),'Py':str(py),'m':str(m),
        'short_tail':str(short),'long_tail':str(long),'c':str(c),'tail_gap':str(1-c),
        'smoothing_halfwidth_max':str(eps),'smoothed_tail_at_max':str(smoothed),
        'smoothed_tail_gap':str(1-smoothed),'guard_x':2,'guard_y':4,
        'gamma_left_break_even_interval':list(map(str,threshold)),
        'gamma_left_gt_3.46882944_is_sufficient_for_positive_K':True,
        'gamma_left_lt_3.46882943_is_sufficient_for_negative_K':True,
        'sequence':'epsilon_j=2^(-j-10), Kx_j=Ky_j=j+9, j>=1; threshold lambda=1 only',
        'source_hashes':{str(source):sha(source),str(cert):sha(cert)},'script_sha256':sha(Path(__file__))}
    dst=base/'computation/results/double-reflection-closure/transfer-independent-check.json'
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','python':sys.version.split()[0],'c_display':float(c),
                     'smoothed_tail_display':float(smoothed),'break_even_display':list(map(float,threshold)),
                     'output':str(dst),'output_sha256':sha(dst)}))
if __name__=='__main__':main()
