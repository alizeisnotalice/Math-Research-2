"""Certify a numerical AH continuation radius using one strict double-hit reference body."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from zipfile import ZipFile

from flint import arb, ctx


def run(comp):
    assert __debug__
    sys.path.insert(0,str(comp))
    from certify_sigma_extended_graph import geometry
    from certify_three_sigma_area_budget import pack, fixed_ah, subtract, exact

    def read(path):
        return json.loads(path.read_bytes())

    def digest(path):
        return sha256(path.read_bytes()).hexdigest()

    baseline_path = comp/'sigma154-turn01150-universal.json'
    baseline = read(baseline_path)
    assert baseline['complete']
    cfg = baseline['parameters']
    pivot, delta0 = map(F,cfg['line'])
    assert (pivot,delta0) == (F('.429'),F('.01150'))
    sigma = cfg['sigma_degrees']
    assert F(sigma)==F('1.54')
    far,negative,eps = map(F,cfg['physical_parameters'])
    ns, src, _ = geometry(comp,[(pivot,delta0)],sigma=sigma,far=far,negative=negative,eps=eps)
    assert sha256(src.encode()).hexdigest() == baseline['chart']['geometry_sha256']
    A,dot,sub,cross,reflect,support = (ns[k] for k in ('A','dot','sub','cross','reflection','supports'))
    w = F('0.05')
    vertex = ns['points'](A(w))[0]
    edges = [sub(vertex[k],vertex[l]) for k,l in ns['orders']]
    # Each cyclic pair is exactly unit by the inverse-circle and DC pairing
    # identities of the construction. Test all remaining edges strictly.
    other_sq = []
    for i,j in combinations(range(5),2):
        u = sub(vertex[i],vertex[j]); ds = dot(u,u)
        if (j-i) in (1,4):
            assert ds.overlaps(A(1))
        else:
            assert ds < A(1)
            other_sq.append(ds.upper())
    assert len(other_sq)==5

    n = list(map(ns['av'],ns['d']['n']))
    m = list(map(ns['av'],ns['d']['m']))
    # L and N lie on DC and DE by their inverse definitions. The second
    # point is the reflected DC mate, so it lies on EF. Only the vertex-fan
    # inclusions are numerical. Opposite supports follow from width one.
    pinned = []
    for j,k in ((3,1),(4,3),(5,0)):
        assert dot(n[j],vertex[k]).overlaps(A(F(1,2)))
        possibilities=[]
        for i,(target,arc) in enumerate(ns['sector']):
            if target != k or arc:
                continue
            lower = min(cross(edges[i],n[j]).lower(),-cross(edges[(i+1)%10],n[j]).upper())
            if lower>0:
                possibilities.append(lower)
        assert len(possibilities)==1
        pinned.append(possibilities[0])
    clock = A(1),A(0)
    cap_room = [A(F(1,2))-support(vertex,edges,clock,u)
                for u in (m[3],m[5],reflect(m[3]),reflect(m[5]))]
    assert min(cap_room)>0
    _, outer, offset = [row for row in ns['rows'] if row[0]=='AH-line:3'][0]
    # These two points are respectively the true straight and curved
    # boundary points of S_A, and hence lie in the same closed base S_A.
    pts = vertex[4], ns['q'](w)
    assert all(u.overlaps(v) for u,v in zip(pts[1],reflect(vertex[2]),strict=True))
    nu = F(1,6000)
    hits = [dot(outer,q)-offset for q in pts]
    assert min(hits) > A(nu)

    archive_path = comp/'sigma154-turn01150-universal.zip'
    ctx.prec = 160
    with ZipFile(archive_path) as z:
        assert z.read('universal-cover.json') == baseline_path.read_bytes()
        selected = {k:v for r in baseline['bins'].values() for k,v in r['sources'].items()}
        assert len(selected)==79
        source_hashes = json.loads(z.read('common-sources.json'))
        for name in ('certify_ah_attainable_root.py','certify_ah_affine_root_class.py',
                     'certify_ah_root_pose_choices.py','certify_ah_pinned_normals.py','certify_sigma_derivative_classes.py'):
            assert digest(comp/name)==source_hashes[name]

        def excluded(node):
            if isinstance(node,dict):
                if node.get('root_free_leaves',0)>0:
                    yield arb(node['sign_separation_lower']).lower(), node['root_free_leaves']
                for child in node.values():
                    yield from excluded(child)
            elif isinstance(node,list):
                for child in node:
                    yield from excluded(child)

        bounds = []
        for name,h in selected.items():
            buf = z.read(name)
            assert sha256(buf).hexdigest()==h
            bounds += [(v,leaves,name) for v,leaves in excluded(json.loads(buf))]
    gap, _, limiting = min(bounds)
    assert gap>0
    root = arb(baseline['root_minimum']).lower()
    pose = arb(baseline['retained_pose_minimum']).lower()

    # The body, its reflection, and the hinge are in the original hexagon.
    # Its diameter is 2/sqrt(3), smaller than rho. The norm of the change
    # in the outer unit normal of AH is at most the change of its angle.
    rho = F(6,5)
    assert exact(rho)>2/arb(3).sqrt()
    increment = F('8e-14')
    t = rho*increment/(nu+rho*increment)
    hausdorff = rho*t
    coeff = F(16)
    root_remaining = root-3*exact(hausdorff)
    gap_remaining = gap-3*exact(hausdorff)
    pose_remaining = pose-exact(coeff*hausdorff+rho*increment)
    assert min(root_remaining,gap_remaining,pose_remaining)>0

    old_ah,crossing0,_ = fixed_ah(comp,sigma,pivot,delta0,far)
    new_ah,crossing1,_ = fixed_ah(comp,sigma,pivot,delta0+increment,far)
    gain = subtract(new_ah,old_ah)
    assert gain[0]>0
    return dict(complete=True, source_sha256=digest(Path(__file__)), baseline_sha256=digest(baseline_path),
                archive_sha256=digest(archive_path),strict_reference_param=str(w),
                strict_reference_fan_margins=[pack(a,24) for a in pinned],
                strict_reference_oblique_margins=[pack(a,24) for a in cap_room],
                strict_reference_nonadjacent_distance_sq_upper=pack(max(other_sq),24),
                strict_reference_AH_hits=[pack(a,24) for a in hits],shared_hit_depth_lower=str(nu),
                original_root_bands_with_exclusions=len(bounds),root_free_leaves=sum(n for _,n,_ in bounds),
                root_exclusion_gap=pack(gap,24),limiting_root_exclusion_source=limiting,
                original_root_margin=pack(root,24),retained_pose_margin=pack(pose,24),
                hausdorff_perturbation_upper=str(hausdorff),normal_term_lipschitz_upper=str(coeff),
                common_AH_angle_band=[str(delta0),str(delta0+increment)],
                root_margin_after_transfer=pack(root_remaining,24),root_exclusion_gap_after_transfer=pack(gap_remaining,24),
                pose_margin_after_transfer=pack(pose_remaining,24),
                deeper_AH_curve_crossing=pack(crossing1,30),round222_AH_curve_crossing=pack(crossing0,30),
                deeper_AH_exact_gain=pack(gain,42),
                newer_area_exact=pack(subtract(tuple(map(F,baseline['eleven_EH']['area']['exact'])),gain),42))


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    p.add_argument('computation',type=Path)
    p.add_argument('output',type=Path)
    a = p.parse_args()
    result = run(a.computation.resolve())
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['common_AH_angle_band'],result['deeper_AH_exact_gain']['decimal'],flush=True)
