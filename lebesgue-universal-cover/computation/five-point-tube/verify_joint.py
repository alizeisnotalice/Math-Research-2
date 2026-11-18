"""Enclose an explicit freely reoptimized seven-input placement with exact rational area bounds."""
import argparse,json,sys
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256

sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_special import area_bound,posed
from template import build


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('comp','points11','pose','output'):p.add_argument(n,type=Path)
    p.add_argument('--k',default=4096,type=int)
    a=p.parse_args();assert __debug__ and a.k>=512
    sys.path.insert(0,str(a.comp.resolve()))
    from certify_nine import append_point,hull
    src=a.points11.read_bytes();assert sha256(src).hexdigest()=='05fb76a261b46c7b939885becb32b31b3e5faca2b017dde5ff70123ff22f6a0d'
    p11=json.loads(src)['points']
    old=(a.comp/'finite-point-leaf-input.json').read_bytes()
    assert sha256(old).hexdigest()=='b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc'
    p30=json.loads(old);p30=[tuple(F(v,p30['denominator']) for v in z) for z in p30['points']]
    data=json.loads(a.pose.read_bytes());pars=list(map(F,data['parameters']))
    assert len(pars)==17 and data['new11_hand'] in (-1,1)
    center=[round(x*2**44) for x in pars[:11]]
    h,D=build(a.comp,False,a.k,True,box=[[x,x] for x in center])
    points=list(h)
    for j,ps,hand in ((11,p30,1),(14,p11,data['new11_hand'])):
        for x,y in posed(ps,dict(pose=[pars[j+1],pars[j+2],pars[j]],hand=hand)):append_point(points,x,y)
    exact=hull(points);lo,hi=area_bound(exact,D,a.k);upper=F((hi*10**12).__ceil__(),10**12)
    result=dict(curve_pose_in_units_2_to_minus44=center,area_lower=str(lo),area_upper=str(hi),upper_decimal=f'{float(upper):.12f}',arc_samples=a.k,
                total_rational_hull_vertices=len(exact),input_sha256=sha256(a.pose.read_bytes()).hexdigest(),new11_sha256=sha256(src).hexdigest(),original_P30_sha256=sha256(old).hexdigest())
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(upper,flush=True)
if __name__=='__main__': main()
