"""Replay a scaled eleven-face EH chain and transfer the proved 1.54-degree AH chart."""
from argparse import ArgumentParser
from contextlib import redirect_stdout
from decimal import Decimal
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys

from flint import arb,ctx

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_three_sigma_area_budget import pack,subtract,enclosure,geometries,setup,ah_clearance
from certify_sigma154_eh_transfer import SIGMA,PIVOT,TURN,D_CHART,exact


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def adapted(f,old,new):
    source=inspect.getsource(f)
    assert source.count(old)==1
    source=source.replace(old,new)
    env=dict(f.__globals__)
    exec(compile(source,f.__name__+'_scaled_for_154','exec'),env)
    return env[f.__name__],sha256(source.encode()).hexdigest()


def run(a):
    comp=a.comp.resolve()
    sys.path.insert(0,str(comp))
    from certify_eh_upper_chain import b
    from certify_eh_upper_chain import front as source_front
    from certify_eh_upper_chain import sign_case as source_sign_case
    original=comp/'eh-two-line-cover-certificate.json'
    params_file=comp/'eh-upper-chain-parameters.json'
    old=load(original)
    params=load(params_file)
    weights=load(a.weights)
    scale=Decimal(weights['scale'])
    assert scale==Decimal('.997')
    adjust=lambda x:format(Decimal(x)*scale,'f')
    old['parameters']=[adjust(d) for d in old['parameters'][:2]]+old['parameters'][2:]
    params=[[adjust(d),p] for d,p in params]
    front,front_hash=adapted(source_front,"F('0.000555')<F(d)","F(old['parameters'][0])<F(d)")
    sign_case,sign_hash=adapted(source_sign_case,"depth='0.000555'","depth=old['parameters'][0]")
    ctx.prec=160
    C,model_sha=b.eh.specialize(F(SIGMA))
    geo=front(old,params)
    g=b.g
    point=b.prev.point_on_faces
    # The first two faces are positive combinations of the two relaxed rows.
    old_norm,offset,r126_gain,u,w=b.eh.clip_geometry(old['parameters'])
    x=point(old_norm,geo['relaxed'])
    t=geo['t']
    assert all(g.cross(g.minus(q,p),g.minus(x,p))>0 for p,q in ((u,w),(w,t),(t,u)))
    det=g.cross(*old_norm)
    combinations=[(g.cross(m,old_norm[1])/det,g.cross(old_norm[0],m)/det)
                  for m in geo['normals'][:2]]
    assert all(a>0 and b>0 for a,b in combinations)
    assert g.cross(g.minus(t,u),g.minus(x,u))<0
    keys=[*old['ordered_branches'],'H0','0H','1H','HH',*(f'U{i}' for i in range(len(params)))]
    assert set(weights['branches'])==set(keys) and len(keys)==17
    bounds={}
    for key in keys:
        p=sign_case(C,old,geo,params,key)
        s=weights['branches'][key]
        b.replay_cuts(p,s['adaptive_cuts'])
        obj=b.eh.phi_objective(p,b.RIGHT,1)
        upper=p.validate(obj,s['dual'])-arb(3)/2
        assert upper<b.arbq(b.F_LIMIT)
        bounds[key]=upper
    depth=params[-1][0]
    g.trial(depth,tight=True)
    b.eh.voronoi_signs(b.eh.arc_data())
    p=C(depth)
    b.eh.add_joint_query(p)
    obj=b.eh.phi_objective(p,b.RIGHT,1)
    p.linear(obj,g.arb(3)/2+b.arbq(b.F_LIMIT))
    replay=weights['pose']
    b.replay_cuts(p,replay['adaptive_cuts'])
    with redirect_stdout(sys.stderr):
        pose=b.eh.certify(p,b.RIGHT,10,replay,eh_relocation=True,joint_eh=True,revisit_ah=True)['summary']
    assert all(v==replay['summary'][k] for k,v in pose.items() if k!='adaptive_cut_rows')
    assert arb(pose['minimum_slack_lower_bound'])>0
    assert arb(pose['F_right_upper'])<0
    full=geometries(SIGMA,old['parameters'],params,g)
    assert full['facets']==geo['active'] and len(full['facets'])==11
    EH_gain=full['new_eh']-full['old_eh']
    alternate=r126_gain+geo['gain']
    assert EH_gain.overlaps(alternate)
    V=setup(SIGMA,g)
    v,t=geo['v'],geo['t']
    arc=g.minus(V['U'],V['C3'])
    c,s=g.unit(exact(D_CHART))
    u555=g.plus(V['C3'],g.minus(g.scale(c,arc),g.scale(s,g.rot(arc))))
    assert v[1]>u555[1]>t[1] and g.cross(g.minus(t,v),g.minus(u555,v))<0
    assert g.cross(g.minus(v,V['C3']),g.minus(u555,V['C3']))>0
    old_n,new_n=b.unit_outer(u555,t),b.unit_outer(v,t)
    cost=g.norm(g.minus(v,t))*g.norm(g.minus(new_n,old_n))
    margin_data=load(a.single)
    assert sha(a.single)=='c3244d638b51cec714d72a4012d807895d262e476c76ca03c52de84419bd3381'
    assert F(margin_data['pivot'])==PIVOT
    assert F(margin_data['EH_chart_depth'])==D_CHART and margin_data['charts']==70
    assert margin_data['sigma']==SIGMA and F(margin_data['AH_turn'])==TURN
    margin=exact(F(margin_data['AH_placement_margin']['exact'][0]))
    assert margin>cost
    clearance=ah_clearance(V,v,t,g,PIVOT,TURN)
    base=tuple(map(F,margin_data['base_area']['exact']))
    ah=tuple(map(F,margin_data['AH_cut']['exact']))
    combined=subtract(subtract(base,ah),enclosure(EH_gain))
    single_area=tuple(map(F,margin_data['combined_area']['exact']))
    reference=comp/'fan-0062-eleven-eh.json'
    published=tuple(map(F,load(reference)['area_exact_interval']))
    saving_single=subtract(single_area,combined)
    saving_published=subtract(published,combined)
    assert saving_single[0]>0 and saving_published[0]>0
    return dict(sigma=SIGMA,pivot=str(PIVOT),AH_turn=str(TURN),scale=str(scale),
                old_parameters=old['parameters'],EH_chain_params=params,EH_deepest=depth,
                EH_11_active_facets=geo['active'],EH_branches=17,
                EH_right_critical_upper=pack(max(bounds.values())),
                EH_F_upper_by_case={k:pack(v) for k,v in bounds.items()},
                room_for_transfer=pack(margin-cost),transfer_cost=pack(cost),
                AH_EH_clearance=pack(clearance),EH_declared_right_bound=pose['F_right_upper'],
                EH_pose_margin=pose['minimum_slack_lower_bound'],EH_cells=pose['angle_edge_cells'],
                EH_duals=pose['validated_duals'],original_baseline=pack(base),AH_cut=pack(ah),
                EH_cut=pack(EH_gain),combined_area=pack(combined),
                saving_from_single=pack(saving_single),saving_from_R201=pack(saving_published),
                source_sha256=sha(Path(__file__)),weights_sha256=sha(a.weights),
                single_sha256=sha(a.single),R126_sha256=sha(original),chain_params_sha256=sha(params_file),
                front_adapter_sha256=front_hash,sign_adapter_sha256=sign_hash,model_sha256=model_sha)


if __name__=='__main__':
    assert __debug__
    parser=ArgumentParser(description=__doc__)
    for key in ('comp','weights','single'):
        parser.add_argument(key,type=Path)
    print(json.dumps(run(parser.parse_args()),separators=(',',':')))
