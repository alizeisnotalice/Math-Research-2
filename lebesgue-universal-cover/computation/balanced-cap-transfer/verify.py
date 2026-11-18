"""Replay all active angle rows and join their shared translation to the projected cap estimate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def certified(comp, work, output):
    here = Path(__file__).resolve().parent
    projected = comp/'projected-cap-transfer'
    rows = module(here/'certify_balanced_rows.py', 'angle_row_certificate')
    cap = module(projected/'verify.py', 'fixed_center_certificate')
    points = module(projected/'generate_joint_cloud.py', 'exact_diameter_clouds')
    work.mkdir(exist_ok=True,parents=True)
    source_cloud, bank, scope = rows.check_paths(comp,work)
    audit = rows.enum(source_cloud,bank,scope)
    assert audit['complete'] and audit['distinct_triples'] == 404
    assert (audit['cells'],audit['row_angle_instances'],audit['represented_handednesses']) == (384,2493,2)
    old = cap.certify(comp,F('.000000872664'))
    assert old == json.loads((projected/'verified.json').read_bytes())
    rate = F(old['projection_rate'])
    q = F(audit['worst_intervals']['fixed_witness'][1])
    eps = 1-F(old['scale'])
    rho = eps*rate/q
    assert F(old['point_shrink']) < rho < 1 and q*rho == eps*rate
    d2 = F(old['point_diameter_squared'])
    factor2 = 1/(d2*(1-rho)**2)
    assert F(old['scale'])**2*factor2 < 1
    bound = factor2*F(old['threshold'])
    with TemporaryDirectory(prefix='balanced-point-replay-',dir=work) as name:
        dest = Path(name)
        points.main(comp,dest,rho,None)
        cloud = (dest/'joint-clouds.json').read_bytes()
        data = json.loads((dest/'joint-receipt.json').read_bytes())
        assert cloud == (here/'joint-clouds.json').read_bytes()
        assert F(data['conditional_area_lower']) == bound and F(data['factor_squared']) == factor2
        assert all(F(row['inner_scale_squared']) >= F(old['scale'])**2*factor2
                   for row in data['group_certificates'][:-1])
        actual = json.loads(cloud)
        checks = []
        for body in actual:
            coords = body['points']
            max_squared = max((x-u)**2+(y-v)**2 for (x,y),(u,v) in combinations(coords,2))
            assert max_squared == body['denominator']**2
            checks.append(len(coords)*(len(coords)-1)//2)
    lower_text = str((bound.numerator*10**24)//bound.denominator).zfill(25)
    decimal = lower_text[:-24]+'.'+lower_text[-24:]
    assert F(decimal) < bound
    result = dict(complete=True,angle_partition=audit,reference_projection_rate=str(rate),
                  point_shrink=str(rho),fixed_point_projection_upper=str(q),
                  transferred_reference_checked_slabs=old['projection_slabs'],
                  transfer_area_lower=str(bound),strict_decimal_lower=decimal,
                  gain_over_231=str(bound-F(old['area_lower'])),direct_pair_checks=sum(checks),
                  all_six_diameter_squared=[1]*6,points_per_input=[len(i['points']) for i in actual],
                  source_sha256={p.name:sha256(p.read_bytes()).hexdigest()
                    for p in (here/'certify_balanced_rows.py',projected/'verify.py',
                              projected/'generate_joint_cloud.py',projected/'verified.json')},
                  point_table_sha256=sha256(cloud).hexdigest())
    output.write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(decimal, result['points_per_input'], sum(checks))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('comp','work','output'):
        p.add_argument(name,type=Path)
    args=p.parse_args()
    assert __debug__
    certified(args.comp.resolve(),args.work.resolve(),args.output.resolve())
