"""Recheck two varying EH faces, leaving the proved AH and deepest EH envelopes fixed."""
from argparse import ArgumentParser
from decimal import Decimal
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, ctx


LAYOUT = {
    4: ('.00000005', '.000000175'),
    7: ('.00000003', '.00000009'),
}


def run(comp, seeds):
    assert __debug__
    sys.path.insert(0,str(comp))
    from certify_eh_upper_chain import b, front as source_front, sign_case as source_sign_case
    from certify_sigma154_eh_chain import adapted
    from certify_three_sigma_area_budget import geometries, setup, pack, enclosure, subtract, exact, ah_clearance

    ctx.prec = 192

    def read(name):
        path = comp/name
        data = path.read_bytes()
        return json.loads(data), sha256(data).hexdigest()

    hull, hull_sha = read('sigma154-turn01150-universal.json')
    source, geometry_sha = read('sigma154-ah429-eh-chain997-composition.json')
    previous, previous_sha = read('eh-two-line-cover-certificate.json')
    weights, weights_sha = read('eh-sigma-154-chain997-weights.json')
    assert geometry_sha == 'd7a2c66a83d64f0afd64ea89ebee6e4ad6c4f1bde292eae286249d00ef8b12ee'
    assert hull['complete'] and hull['parameters']['line'] == ['.429','.01150']
    assert weights_sha == source['weights_sha256'] and previous_sha == source['R126_sha256']
    assert hull['parameters']['sigma_degrees'] == source['sigma'] == '1.54'
    assert F(hull['eleven_EH']['EH_deepest']) == F(source['EH_deepest'])
    pivot, turn = map(F,hull['parameters']['line'])
    old = {**previous, 'parameters':source['old_parameters']}
    base = source['EH_chain_params']
    g = b.g
    model, model_sha = b.eh.specialize(F('1.54'))
    assert source['model_sha256'] == model_sha
    front, front_sha = adapted(source_front,"F('0.000555')<F(d)","F(old['parameters'][0])<F(d)")
    sign_case, sign_sha = adapted(source_sign_case,"depth='0.000555'","depth=old['parameters'][0]")
    assert front_sha == source['front_adapter_sha256'] and sign_sha == source['sign_adapter_sha256']
    og = geometries('1.54',old['parameters'],base,g)
    reference = front(old,base)
    setup_coords = setup('1.54',g)
    n, _, _, u, _ = b.eh.clip_geometry(old['parameters'])
    x = b.prev.point_on_faces(n,reference['relaxed'])
    T,C3,U = [setup_coords[key] for key in ('T','C3','U')]
    radial = g.minus(U,C3)

    def facet_pivot(j):
        return g.plus(x,g.scale(exact(base[j][1]),g.minus(u,x)))

    def corner(d):
        c,s = g.unit(d)
        return g.plus(C3,g.minus(g.scale(c,radial),g.scale(s,g.rot(radial))))

    def decimal(d):
        r = format(Decimal(d.numerator)/Decimal(d.denominator),'f')
        assert F(r)==d
        return r

    def parameter_set(increments):
        params = [list(pair) for pair in base]
        for j,s in increments.items():
            params[j][0] = decimal(F(params[j][0])+F(s))
        return params

    upper = {j:pair[1] for j,pair in LAYOUT.items()}
    lower = {j:pair[0] for j,pair in LAYOUT.items()}
    low_params,high_params = parameter_set(lower),parameter_set(upper)
    lg,hg = (geometries('1.54',old['parameters'],params,g) for params in (low_params,high_params))
    f = front(old,high_params)
    assert og['facets'] == lg['facets'] == hg['facets'] == f['active'] == source['EH_11_active_facets']
    assert len(hg['facets']) == 11
    V = corner(exact(base[-1][0]))
    clearance = ah_clearance(setup_coords,V,T,g,pivot,turn)
    assert clearance > 0
    # The baseline and all perturbed envelopes use the same upper endpoint.
    assert base[-1] == low_params[-1] == high_params[-1]

    branch_upper = {}
    slope_increase = {}
    seed_hashes = {}
    for j,seed_path in zip(LAYOUT,seeds,strict=True):
        seed = json.loads(seed_path.read_bytes())
        key = f'U{j}'
        assert seed['facet_index']==j+2 and seed['changed_branch']==key
        assert seed['base_EH_geometry_sha256']==geometry_sha and seed['base_weights_sha256']==weights_sha
        selection, = [r for r in seed['rows'] if r['delta']==upper[j] and r['passed']]
        data = selection['branch_witness']
        p = sign_case(model,old,f,high_params,key)
        b.replay_cuts(p,data['adaptive_cuts'])
        upper_f = p.validate(b.eh.phi_objective(p,b.RIGHT,1),data['dual'])-arb(3)/2
        assert upper_f < exact(b.F_LIMIT)
        branch_upper[key] = pack(upper_f,24)
        seed_hashes[seed_path.name] = sha256(seed_path.read_bytes()).hexdigest()
        y = facet_pivot(j)
        d0,d1 = F(base[j][0]),F(high_params[j][0])
        ray = corner(exact(d0).union(exact(d1)))
        v = g.minus(ray,y)
        assert v[1] > 0
        dv = g.neg(g.rot(g.minus(ray,C3)))
        slope = g.cross(dv,v)/v[1]**2
        assert slope > 0
        for d in (d0,d1):
            normal = b.unit_outer(corner(exact(d)), y)
            assert normal[0] < 0 and g.cross(g.minus(y,corner(exact(d))),g.minus(x,y)) > 0
        slope_increase[key] = pack(slope,24)

    # Choose one interior witness on each true facet. For a moved facet the
    # hardest corner holds that face shallow and all its competitors deep.
    witnesses = {}
    triangle = [V,U,T]
    for label in hg['facets']:
        j = label-2
        if j in LAYOUT:
            adverse = {**upper,j:lower[j]}
            shell = geometries('1.54',old['parameters'],parameter_set(adverse),g)
        else:
            shell = hg
        assert shell['facets'] == hg['facets']
        k = shell['facets'].index(label)
        a,bpoint = shell['points'][k:k+2]
        z = g.scale(arb(1)/2,g.plus(a,bpoint))
        if j in LAYOUT:
            y = facet_pivot(j)
            move = corner(exact(high_params[j][0]))
            width = (move[0]-y[0])/(move[1]-y[1])
            assert z[1] > y[1]
            last = (y[0]+width*(z[1]-y[1]),z[1])
            assert all(g.cross(g.minus(r,p),g.minus(last,p))>0 for p,r in zip(triangle,triangle[1:]+triangle[:1]))
        witnesses[str(label)] = pack(z[1],24)

    old_area = tuple(map(F,hull['eleven_EH']['area']['exact']))
    area_gain = {}
    area = {}
    for name,geo in (('least_cut_corner',lg),('most_cut_corner',hg)):
        gain = enclosure(geo['new_eh']-og['new_eh'])
        assert gain[0] > 0
        area_gain[name] = pack(gain,42)
        area[name] = pack(subtract(old_area,gain),42)
    return dict(complete=True,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                baseline_round222_sha256=hull_sha, EH_source_sha256=geometry_sha, EH_first_pose_weights_sha256=weights_sha,
                changed_branch_sources=seed_hashes,baseline_EH_depths=base,low_corner=low_params,high_corner=high_params,
                first_stage_changed_branch_upper=branch_upper,continuous_moved_support_slope_lower=slope_increase,
                active_facets=hg['facets'],facet_witness_heights=witnesses,deepest_envelope=source['EH_deepest'],
                EH_AH_separation=pack(clearance,25),area_gain=area_gain,area=area)


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('computation',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('seeds',type=Path,nargs=2)
    a=parser.parse_args()
    result=run(a.computation.resolve(),a.seeds)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result['area'].items()},indent=2))
