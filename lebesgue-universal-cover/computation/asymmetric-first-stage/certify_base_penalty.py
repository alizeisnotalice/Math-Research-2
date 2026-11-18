"""Area cost of unequal caps compared with the larger of their two matched-angle base hosts."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
from flint import arb,ctx
def run(result):
    ctx.prec=192;delta=arb.pi()*arb('.04')/180
    s=arb('1.50').union(arb('1.54'))*arb.pi()/180;r=arb(3).sqrt();z=s.cos()
    match_der=2*s.sin()*(2*r*z-3)**2/(3*(4*z**2-3)**2)
    k=min(arb(row['coefficient_lower']).lower() for row in json.loads(result.read_bytes())['cases'])
    assert k>arb('.11858')
    mass=arb('.11858')*delta
    d_area=2*r*(z/r-arb('.5')-mass)/(4*z*z-3)
    vs_mean=d_area.lower()*arb('.11858')
    vs_larger=vs_mean-match_der.upper()/2
    assert vs_larger>arb('.0298')
    sig3=arb('1.50')*arb.pi()/180;sig5=arb('1.54')*arb.pi()/180
    t=(arb.pi()/12-sig3/2).tan()/2
    add=t*delta.sin()-(1-delta.cos())/2
    def cap(s,e):
        c=s.cos()
        return r*(c/r-arb('.5')-e)**2/(4*c*c-3)
    upper=r/2-cap(sig3,arb(0))-cap(sig5,add)
    assert upper<arb('.845373720')
    return dict(complete=True,specific_sharp_inherited_upper=str(upper.upper()),necessary_cap_height_sum_coefficient='.11858',max_height_for_derivative=str(mass.upper()),
        cap_area_derivative_lower=str(d_area.lower()),matched_area_derivative_upper=str(match_der.upper()),
        penalty_vs_average_per_radian=str(vs_mean.lower()),penalty_vs_larger_per_radian=str(vs_larger.lower()),
        dependency=sha256(result.read_bytes()).hexdigest(),verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser();p.add_argument('data',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();r=run(a.data);a.out.write_text(json.dumps(r,indent=1)+'\n');print(r['penalty_vs_larger_per_radian'])
