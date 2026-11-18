"""A mandatory large support excursion at 65/128 for every host inside P retaining the .95 deep line."""
import argparse,inspect,json,sys
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import certify_ah_released_tangent_obstruction as obstruction
from flint import arb


def run(archive):
    old=inspect.getsource(obstruction.run)
    before="A('.0002') if long_patch"
    t=F('.00335462')
    assert old.count(before)==1
    tested=old.replace(before,f"A('{t}') if long_patch")
    scope=dict(vars(obstruction))
    exec(compile(tested,'explicit_larger_excursion','exec'),scope)
    j=scope['run'](archive,True)
    # run() checks the identity, ten genuinely excluded P caps, and all non-dihedral angles;
    # its QA residual is against the OLD long chain's 1e-6 offset.
    u=arb(j['QA_violation'])+arb('.000001')
    assert u>arb('.0000289') and arb(j['identity_violation'])>0
    return dict(w='65/128',weight=str(t),deep_pivot='.950',deep_angle='.011',
                necessary_outward_displacement_from_old_q_tangent=str(u),
                identity_deep_violation=j['identity_violation'],ten_other_caps=j['ten_other_pose_min_violations'],
                angular_cells=j['interior_angular_cells'],interior_negative_margin=j['interior_negative_margin'],
                endpoint_slope_margin=j['endpoint_slope_clearance'],specialized_source_sha256=sha256(tested.encode()).hexdigest(),
                base_sha256=sha256(Path(obstruction.__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('output',type=Path)
    opt=cli.parse_args()
    path=opt.archive.resolve()
    sys.path.insert(0,str(path))
    result=run(path)
    opt.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
