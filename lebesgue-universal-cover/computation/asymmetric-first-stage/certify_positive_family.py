"""Exact simultaneous offset budget and hull areas for an asymmetric pair of caps."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from flint import arb,ctx
def run():
    ctx.prec=192
    rad=arb.pi()/180
    b=(arb('14.25')*rad).tan()/2
    def displacement(d):
        d=arb(str(d))*rad
        return b*d.sin()+arb(1)/2*(1-d.cos())
    def area(sig3,sig5,e3,e5):
        # Isolated corner at radius 1/sqrt(3), with cap normal tilt sigma.
        rt=arb(3).sqrt()
        def cut(sig,offs):
            angle=arb(str(sig))*rad;c=angle.cos()
            return rt*(c/rt-(arb(1)/2+offs))**2/(4*c**2-3)
        return rt/2-cut(sig3,e3)-cut(sig5,e5)
    e=displacement('.02');e1=displacement('.04')
    tests={}
    for k,e3,e5 in [('equal height',e,e),('only second released',arb(0),e1),('only first released',e1,arb(0))]:
        val=area('1.50','1.54',e3,e5)
        tests[k]=dict(e3=str(e3),e5=str(e5),area=str(val))
    ideal=area('1.50','1.54',arb(0),arb(0));current=arb('.844089981544184759234377880873755496061514')
    assert ideal>current and min(area('1.50','1.54',x,y).lower() for x,y in [(e,e),(arb(0),e1),(e1,arb(0))])>ideal
    return dict(complete=True,family=['1.50','1.54'],radial_width_coefficient=str(b),unreleased_area=str(ideal),
        examples=tests,verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args();r=run()
    a.out.write_text(json.dumps(r,indent=1)+'\n');print(r['examples'],r['unreleased_area'])
