"""Certify full-body containment in the Gibbs hexagon and original two caps.

The six hexagon contacts are exact for the constructions. Arb checks the
pentagon sectors identifying those contacts and strict signed support gaps at
the caps. The script also certifies witnesses against prior walls.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

from flint import __version__ as flint_version, arb, arb_series, ctx


SIGMA = '1.5497'
PENTAGONS = ('nonideal-pentagon-library.json', 'nonmonotone-pentagon-inputs.json')
PROFILES = ('tapered-forced-support-profiles.json', 'nonmonotone-curvature-profiles.json',
            'tapered-double-deep-profile.json')
MIRROR = 'ah-tapered-mirror-counterexample-input.json'
COUNTS = dict(zip(PENTAGONS + PROFILES + (MIRROR,), (105, 12, 19, 11, 1, 2), strict=True))
CONVENTION = ('row vectors: (x,y) @ [[cos(k),s*sin(k)],[-sin(k),s*cos(k)]], '
              'k in degrees; source indices zero-based')


def exact(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def unit(angle):
    sn, co = angle.sin_cos()
    return co, sn


def preimage(normal, turn, sign):
    co, sn = unit(turn * arb.pi() / 180)
    x, y = normal
    return co*x + sn*sign*y, -sn*x + co*sign*y


def transform(point, turn, sign):
    co, sn = unit(turn * arb.pi() / 180)
    x, y = point
    return co*x - sn*y, sign*(sn*x + co*y)


def hull_of_values(values):
    values = iter(values)
    result = next(values)
    for value in values:
        result = result.union(value)
    return result


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_pentagon(ns, vertices, support, hex_normals):
    orders, sector = ns['orders'], ns['sector']
    edges = [ns['sub'](vertices[i], vertices[j]) for i, j in orders]
    for k, (edge, other) in enumerate(zip(edges, edges[1:] + edges[:1], strict=True)):
        assert edge[1] > 0 if k < 5 else edge[1] < 0
        assert ns['cross'](edge, other) > 0
        assert ns['dot'](edge, other) > 0
        i, shift = sector[k]
        assert orders[k][shift] == orders[(k+1) % 10][shift] == i

    def sectors(u):
        tests = [ns['cross'](edge, u) for edge in edges]
        return {sector[k] for k, value in enumerate(tests)
                if not value < 0 and not tests[(k+1) % 10] > 0}

    # Exact pinned contacts follow from L on n5, N on n4, and Q_A L' on n3.
    # Their antipodal contacts follow from width one.
    for j, normal in enumerate(hex_normals):
        assert sectors(normal) == {((1, 3, 0)[j % 3], int(j < 3))}
        assert support(normal).contains(exact('1/2'))

    def support_point(u):
        points = [ns['add'](vertices[i], ns['mul'](shift, u)) for i, shift in sorted(sectors(u))]
        assert points
        return tuple(hull_of_values(coords) for coords in zip(*points, strict=True))

    return support_point


def mirror_support_point(mirror, u, base, g):
    pi = arb.pi()
    relative = arb.atan2(u[1], u[0]) - pi/6
    if relative < -pi:
        relative += 2*pi
    assert not relative.contains(0)
    relative_sign = 1 if relative > 0 else -1
    theta = abs(relative)
    width_sign = 1
    if theta > pi/2:
        theta = pi-theta
        width_sign = -1
    assert 0 < theta < pi/2
    sin_theta, cos_theta = theta.sin_cos()
    breaks = base.angles if mirror is base else mirror.breaks
    jets = [(exact('1/2') + width_sign*(a*cos_theta+b*sin_theta+d),
             relative_sign*(b*cos_theta-a*sin_theta))
            for lo, hi, (a, b, d) in zip(breaks[:-1], breaks[1:], mirror.rows, strict=True)
            if not (theta < lo or theta > hi)]
    assert jets
    h, derivative = (hull_of_values(part) for part in zip(*jets, strict=True))
    return g.plus(g.scale(h, u), g.scale(derivative, g.rot(u)))


def run(computation, poses):
    # Geometry loaders inherited from the project set their own contexts.
    sys.path.insert(0, str(computation))
    from certify_ah_tapered_minkowski_closure import geometry
    from certify_ah_tapered_mirror_counterexample import Polygon, Support, compatible_angles
    from certify_tapered_forced_support import Profile
    from deformed_reuleaux_support import body, environment

    old, old_geometry_hash = geometry(computation)
    ns = environment(str(computation), SIGMA)
    ctx.prec = 192
    g, pi = old['g'], arb.pi()
    half = exact('1/2')
    hex_normals = [unit(j*pi/3) for j in range(6)]
    caps = {j: unit((30 + exact(SIGMA) + 60*j)*pi/180) for j in (3, 5)}
    for j, u in caps.items():
        assert all(x.overlaps(y) for x, y in zip(u, old['tilt'][j], strict=True))

    raw = json.loads(poses.read_bytes())
    assert raw['convention'] == CONVENTION
    records = raw['poses']
    assert sorted(r['index'] for r in records) == list(range(150))
    assert Counter(r['source'] for r in records) == COUNTS
    for source, count in COUNTS.items():
        assert sorted(row['body_index'] for row in records if row['source'] == source) == list(range(count))
    inputs = {name: json.loads((computation / name).read_bytes()) for name in COUNTS}
    mirror = Support(inputs[MIRROR], pi)
    polygon = Polygon(mirror)
    mirror_objs = (mirror, polygon)

    direct, z, D = (old['ns'][key] for key in ('direct', 'z', 'D'))
    original_ah_normal = g.neg(g.rot(direct))
    original_ah = (original_ah_normal, g.dot(original_ah_normal, z))
    angle = exact('.006')
    u = g.neg(g.rot(g.plus(g.scale(angle.cos(), direct), g.scale(angle.sin(), g.rot(direct)))))
    r199_ah = (u, g.dot(u, g.plus(z, g.scale(exact('.427'), D))))
    q = old['q'](arb_series([exact('59/128'), 1]))
    point, velocity = tuple(component[0] for component in q), tuple(component[1] for component in q)
    normal = g.scale(-1/g.norm(velocity), g.rot(velocity))
    paired = (normal, g.dot(normal, point))
    r203 = {name:(u, h) for name, u, h in old['active']}
    r203_taper = r203['taper_23']
    walls = {'original-Ta-Z': original_ah, 'R199-AH': r199_ah,
             'original-paired-59/128': paired, 'R203-taper-23': r203_taper,
             'R203-deep-AH': r203['deep']}
    assert g.dot(paired[0], z) < paired[1]
    required_witnesses = {113: ('original-paired-59/128', 'R203-taper-23'),
                          117: ('R203-deep-AH',), 118: ('R203-deep-AH',)}
    compatibility = compatible_angles(mirror)
    assert 3*polygon.error < exact('.003')
    mirror_exclusion = {}
    per_body, witnesses = [], []
    for row in records:
        source, index = row['source'], row['body_index']
        turn, sign = row['turn'], row['reflection']
        assert type(turn) is int and turn in range(0, 360, 60) and sign in (-1, 1)
        assert all(any(all(a.overlaps(b) for a,b in zip(n, preimage(u, turn, sign), strict=True))
                       for n in hex_normals) for u in hex_normals)
        if source in PENTAGONS:
            if source == PENTAGONS[0]:
                inp = inputs[source]['proofs'][index]
                w, e = inp['w'], inp['deformation']
            else:
                w, e = inputs[source][index]
            w, e = Fraction(w), Fraction(e)
            assert w == Fraction(row['w']) and e == Fraction(row['e'])
            vertices, support = body(ns, w, e)
            support_point = check_pentagon(ns, vertices, support, hex_normals)
            hex_proof = 'exact contact lines and verified ten-sector fan'
        elif source == MIRROR:
            obj = mirror_objs[index]
            assert row['body_type'] == type(obj).__name__
            support = obj.support_direction
            def support_point(u):
                return mirror_support_point(obj, u, mirror, g)
            hex_proof = 'exact mirror symmetry, zero moment and centered matching moments'
            exclusions = []
            for k, face in inputs[MIRROR]['rotation_faces']:
                u, h = r203[face]
                gap = support(preimage(u, k, 1))-h
                assert gap > exact('.00001901')
                exclusions.append({'turn': k, 'face': face, 'excess': str(gap)})
            assert [record['turn'] for record in exclusions] == list(range(0, 360, 60))
            mirror_exclusion[row['index']] = exclusions
        else:
            obj = Profile(inputs[source][index])
            support = obj.support
            hex_proof = 'exact defining linear system for the normalized support'

        def posed_support(normal):
            return support(preimage(normal, turn, sign))

        margins = {str(k): half-posed_support(u) for k,u in caps.items()}
        assert all(gap > 0 for gap in margins.values()), (row['index'], margins)
        ah_margin = original_ah[1]-posed_support(original_ah[0])
        assert ah_margin > 0
        r199_margin = r199_ah[1]-posed_support(r199_ah[0])
        per_body.append((row['index'], margins, ah_margin, r199_margin, hex_proof))
        for label in required_witnesses.get(row['index'], ()):
            u, h = walls[label]
            pulled = preimage(u, turn, sign)
            x = transform(support_point(pulled), turn, sign)
            excess = g.dot(u,x)-h
            assert excess > 0 and excess.overlaps(posed_support(u)-h)
            witnesses.append({'index': row['index'], 'wall': label, 'point': list(map(str, x)),
                              'normal': list(map(str, u)), 'offset': str(h), 'excess': str(excess)})

    def minimum(getter):
        index, value = min(((t[0], getter(t)) for t in per_body), key=lambda q:q[1].lower())
        return {'index': index, 'margin': str(value)}

    deps = tuple(COUNTS) + ('certify_ah_tapered_minkowski_closure.py', 'certify_ah_tapered_release.py',
                            'ah-tapered-release.json', 'certify_ah_tapered_mirror_counterexample.py',
                            'certify_tapered_forced_support.py', 'deformed_reuleaux_support.py',
                            'certify_sigma_extended_graph.py', 'certify_ah_two_line_cap.py',
                            'certify_gibbs_area_refined.py')
    return {'complete': True, 'precision_bits': ctx.prec, 'python': platform.python_version(),
            'python_flint': flint_version, 'sigma_degrees': SIGMA, 'body_count': len(records),
            'all_bodies_exactly_touch_all_six_hex_faces': True,
            'pose_counts': [[k,s,c] for (k,s),c in sorted(Counter((r['turn'],r['reflection']) for r in records).items())],
            'minimum_cap_margins': {str(k): minimum(lambda t: t[1][str(k)]) for k in caps},
            'minimum_original_Ta_Z_margin': minimum(lambda t: t[2]),
            'minimum_R199_AH_margin': minimum(lambda t: t[3]),
            'per_body': [{'index': i, 'cap_margins': {k:str(v) for k,v in margins.items()},
                          'original_Ta_Z_margin': str(a), 'R199_AH_margin': str(b),
                          'hex_proof': proof} for i, margins, a, b, proof in per_body],
            'crossing_points': witnesses,
            'R207_counterexamples_excluded_by_R203': {'mirror_reduced_angles': compatibility,
                                                      'polygon_threefold_error_bound': str(3*polygon.error),
                                                      'old_wall_exclusions': mirror_exclusion},
            'poses_sha256': digest(poses), 'derived_old_geometry_sha256': old_geometry_hash,
            'dependencies_sha256': {k:digest(computation/k) for k in deps},
            'script_sha256': digest(Path(__file__))}


def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    parser = ArgumentParser(description=__doc__)
    for name in ('computation', 'poses', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    result = run(args.computation.resolve(), args.poses.resolve())
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS', result['body_count'], result['minimum_cap_margins'],
          result['minimum_original_Ta_Z_margin'], result['minimum_R199_AH_margin'])
    for row in result['crossing_points']:
        print(row['index'], row['wall'], row['excess'])


if __name__ == '__main__':
    main()
