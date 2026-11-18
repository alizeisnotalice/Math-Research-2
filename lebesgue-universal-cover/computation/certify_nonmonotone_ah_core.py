"""Lower bound the area of universal convex subsets of P containing the retained core."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import __version__ as flint_version, arb, ctx

from certify_ah_tapered_minkowski_closure import geometry
from certify_ah_tapered_mirror_counterexample import Polygon, Support
from certify_nonideal_pentagon_pose_library import verify
from certify_tapered_forced_support import Profile, rational
from deformed_reuleaux_support import body, environment


SCALE = 10**11
ANCHORS = ((0, 0), (SCALE, 0))
ROOT_POSES = [(k, s) for k in range(0, 360, 60) for s in (1, -1)]
BASE_LO = F('.8439741301130277')
LENGTH2_LO = F('.01007656807905897')
INPUTS = (
    'nonideal-pentagon-library.json', 'nonmonotone-pentagon-inputs.json',
    'ah-tapered-mirror-counterexample-input.json', 'ah-tapered-mirror-counterexample-audit.json',
    'tapered-forced-support-profiles.json', 'nonmonotone-curvature-profiles.json',
    'tapered-double-deep-profile.json', 'ah-flexible-chord-certificate.json',
    'flexible-006-eleven-eh.json', 'certify_ah_tapered_minkowski_closure.py',
    'certify_ah_tapered_mirror_counterexample.py', 'certify_tapered_forced_support.py',
    'certify_nonideal_pentagon_pose_library.py', 'deformed_reuleaux_support.py',
)


def hull(points):
    def half(ordered):
        edge = []
        for x, y in ordered:
            while len(edge) > 1:
                a, b = edge[-2:]
                if (b[0]-a[0])*(y-b[1]) > (b[1]-a[1])*(x-b[0]):
                    break
                edge.pop()
            edge.append((x, y))
        return edge

    unique = sorted(set(points))
    if len(unique) < 3:
        return tuple(unique)
    return tuple(half(unique)[:-1] + half(reversed(unique))[:-1])


def double_area(poly):
    return sum(x*q - y*p for (x,y),(p,q) in zip(poly, poly[1:] + poly[:1]))


def profile_jet(profile, normal):
    theta = arb.atan2(normal[1], normal[0])
    pi = arb.pi()
    if theta < 0:
        theta += 2*pi
    sign = 1
    if theta > pi:
        theta -= pi
        sign = -1
    assert 0 < theta < pi
    for k, curvature in enumerate(profile.fs):
        if not profile.knots[k] < theta < profile.knots[k+1]:
            continue
        sn, co = (theta - profile.knots[k]).sin_cos()
        y, v = profile.states[k]
        return arb(1)/2 + sign*(curvature+(y-curvature)*co+v*sn), sign*(v*co-(y-curvature)*sn)
    return None


def mirror_jet(mirror, normal):
    pi = arb.pi()
    theta = arb.atan2(normal[1], normal[0])-pi/6
    if theta < -pi:
        theta += 2*pi
    assert not theta.contains(0)
    mirror_sign = 1 if theta > 0 else -1
    q = abs(theta)
    width_sign = 1
    if q > pi/2:
        q = pi-q
        width_sign = -1
    assert 0 < q < pi/2
    breaks = mirror.angles if isinstance(mirror, Support) else mirror.breaks
    for k, (l, r) in enumerate(zip(breaks, breaks[1:])):
        if l < q < r:
            a,b,d = mirror.rows[k]
            sn,co = q.sin_cos()
            return arb(1)/2+width_sign*(a*co+b*sn+d), mirror_sign*(b*co-a*sn)
    return None


def prepare(comp, steps, cap_shift):
    old,_ = geometry(comp)
    g, source = old['g'], old['ns']
    ns = environment(str(comp), '1.5497')
    ctx.prec = 160
    ta = g.plus(source['z'],source['D'])
    tp = old['q'](arb(1))
    v = g.minus(tp,ta)
    v2 = g.dot(v,v)
    core = old['oldarea']+old['oldcap']-old['green'](F(0),F(1),ta)
    assert core > rational(BASE_LO) and v2 > rational(LENGTH2_LO)
    assert g.cross(v,g.neg(ta)) > 0 and g.cross(v,g.minus(source['z'],ta)) < 0
    outside = g.neg(g.rot(v))
    normals = []
    for i in range(-steps,steps+1):
        n = g.plus(outside, g.scale(arb(i)/(steps*3),v))
        normals.append(g.scale(1/g.norm(n),n))

    def inner_point(point):
        delta = g.minus(point,ta)
        x = g.dot(v,delta)/v2
        y = -g.cross(v,delta)/v2
        if not y > arb(10)/SCALE:
            return None
        # A rational point at this height lies between the two triangle sides.
        iy = int(float(y.mid())*SCALE)-4
        t = arb(iy)/(SCALE*y)
        assert 0 < t < 1
        left,right = t*x, t*x+1-t
        ix = round(float(((left+right)*SCALE/2).mid()))
        assert left*SCALE < ix < right*SCALE
        return ix,iy

    def cap_options(poses, jet):
        opts = set()
        for k,s in poses:
            points = set(ANCHORS)
            for n in normals:
                pre = old['preimage'](n,k,s)
                value = jet(pre)
                if value is None:
                    continue
                h, derivative = value
                # det(T)=s changes the sign of the derivative in the host frame.
                point = g.plus(g.scale(h,n),g.scale(s*derivative,g.rot(n)))
                p = inner_point(point)
                if p is not None:
                    points.add(p)
            opts.add(hull(points))
        assert opts
        return sorted(opts)

    caps = [old['tilt'][i] for i in (3,5)]
    cap_margins = []
    def allowed(h, poses):
        kept = []
        for k,s in poses:
            gap = max((h(old['preimage'](n,k,s))-arb(1)/2).lower() for n in caps)
            if gap > rational(cap_shift):
                cap_margins.append(gap)
            else:
                # An indeterminate sign retains the pose in this relaxation.
                kept.append((k,s))
        return kept

    groups = []
    def load(name):
        return json.loads((comp/name).read_bytes())
    library = load('nonideal-pentagon-library.json')
    assert not library['unresolved'] and len(library['proofs']) == 105
    extra = [verify(ns,F(w),F(e),F(0),F(0)) for w,e in load('nonmonotone-pentagon-inputs.json')]
    for index,row in enumerate(library['proofs']+extra):
        vertices,support = body(ns,F(row['w']),F(row['deformation']))
        edges = [g.minus(vertices[i],vertices[j]) for i,j in ns['orders']]
        assert {(z[0],z[1]) for z in row['cap_eliminations']} == set(ROOT_POSES)-{(0,1),(240,-1)}
        assert all(arb(z[3])>0 for z in row['cap_eliminations'])
        poss = [(0,1),(240,-1)]
        for k,s,_,saved_gap in row['cap_eliminations']:
            if arb(saved_gap) > rational(cap_shift):
                cap_margins.append(arb(saved_gap).lower())
            else:
                poss.extend(allowed(support,[(k,s)]))
        assert arb(row['interior_margin'])>0 and arb(row['endpoint_margin'])>0
        def pentagon_jet(n):
            for j,e in enumerate(edges):
                if g.cross(e,n)>0 and g.cross(n,edges[(j+1)%len(edges)])>0:
                    i,r = ns['sector'][j]
                    return g.dot(vertices[i],n)+r, g.dot(vertices[i],g.rot(n))
            return None
        groups.append((f'pentagon-{index}',cap_options(poss,pentagon_jet)))

    mirror = Support(load('ah-tapered-mirror-counterexample-input.json'),arb.pi())
    mirror_proof = load('ah-tapered-mirror-counterexample-audit.json')
    assert mirror_proof['complete'] and mirror_proof['input_sha256'] == sha256((comp/'ah-tapered-mirror-counterexample-input.json').read_bytes()).hexdigest()
    for label, curve in [('mirror',mirror),('mirror-polygon',Polygon(mirror))]:
        # Mirror symmetry was proved with its 30-degree axis in the referenced proof.
        poss = allowed(curve.support_direction,[(k,1) for k in range(0,360,60)])
        groups.append((label,cap_options(poss,lambda n:mirror_jet(curve,n))))

    profile_sources = ['tapered-forced-support-profiles.json','nonmonotone-curvature-profiles.json','tapered-double-deep-profile.json']
    for name in profile_sources:
        for i,spec in enumerate(load(name)):
            curve = Profile(spec)
            curve.root_clearance()
            poss = allowed(curve.support,ROOT_POSES)
            groups.append((f'{name}/{i}',cap_options(poss,lambda n:profile_jet(curve,n))))
    margin = min(cap_margins)
    # The interval comparison, after selecting a decimal, certifies this common relaxation.
    eta = F(int(float(margin.mid())*10**15)-1,10**15)
    eta = max(eta,cap_shift)
    assert 0 < eta and margin > rational(eta)
    return groups,{'core_area':str(core),'chord_squared':str(v2),'normal_count':len(normals),
                   'original_cap_gap_lower':str(margin),'allowed_outward_shift':str(eta)}


def descend(groups, target):
    threshold = (target-BASE_LO)*2*SCALE**2/LENGTH2_LO
    groups = sorted(groups,key=lambda t:min(double_area(h) for h in t[1]),reverse=True)
    states = {ANCHORS}
    layers = []
    min_pruned = None
    for name,options in groups:
        next_states = set()
        for state in states:
            enlarged = {hull(state+choice) for choice in options}
            if state in enlarged:
                enlarged = {state}
            for outcome in enlarged:
                a = double_area(outcome)
                if a*threshold.denominator > threshold.numerator:
                    min_pruned = a if min_pruned is None else min(a,min_pruned)
                else:
                    next_states.add(outcome)
        layers.append([name,len(options),len(states),len(next_states)])
        states = next_states
        if not states or len(states)>4096:
            break
    least = min(([min_pruned] if min_pruned is not None else [])+[double_area(s) for s in states])
    lower = BASE_LO+LENGTH2_LO*least/(2*SCALE**2)
    return {'target':str(target),'guaranteed_lower':str(lower),'states_below_target':len(states),'layers':layers}


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('computation',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--target',type=F,default=F('.844081861'))
    parser.add_argument('--steps',type=int,default=540)
    parser.add_argument('--cap-shift',type=F,default=F(0))
    args = parser.parse_args()
    assert args.steps > 0 and args.cap_shift >= 0
    groups,rec = prepare(args.computation.resolve(),args.steps,args.cap_shift)
    result = descend(groups,args.target)
    rec.update(result,group_count=len(groups),grid_scale=SCALE,steps=args.steps,flint=flint_version,
               requested_cap_shift=str(args.cap_shift),
               input_sha256={n:sha256((args.computation/n).read_bytes()).hexdigest() for n in INPUTS},
               source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(rec,indent=2)+'\n')
    print(len(groups),len(result['layers']),result['states_below_target'],result['guaranteed_lower'],flush=True)
