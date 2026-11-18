"""Replay the first EH placement with a scaled eleven-face chain at a given angle."""
from argparse import ArgumentParser
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from certify_sigma154_eh_chain import adapted
from certify_three_sigma_area_budget import enclosure, geometries, pack

BASE_GEOMETRY = (
    "certify_ah_shortest.py", "certify_ah_contact_order.py", "audit_eh_geometry_final.py",
    "certify_splice_connectors.py", "certify_cover_relocations.py",
)


def run(comp, sigma, weights):
    sys.path.insert(0,str(comp))
    from certify_eh_upper_chain import b, front as base_front, sign_case as base_sign

    original_env = {**os.environ, "LEBESGUE_SIGMA_DEGREES": f"{sigma},{sigma}"}
    checks = {}
    for name in BASE_GEOMETRY:
        raw = subprocess.check_output([sys.executable,"-E","-B",str(comp/name)], env=original_env)
        d = json.loads(raw)
        assert list(map(F,d["sigma_degrees_interval"])) == [F(sigma),F(sigma)]
        assert d.get("passed",d.get("exact_maps_passed")) is True
        checks[name] = sha256(raw).hexdigest()

    frozen = json.loads((HERE/"sigma154-ah429-eh-chain997-composition.json").read_bytes())
    model_old = json.loads((comp/"eh-two-line-cover-certificate.json").read_bytes())
    model_old["parameters"] = frozen["old_parameters"]
    params = frozen["EH_chain_params"]
    front, front_hash = adapted(base_front,"F('0.000555')<F(d)","F(old['parameters'][0])<F(d)")
    sign_case, sign_hash = adapted(base_sign,"depth='0.000555'","depth=old['parameters'][0]")
    ctx.prec=160
    C, model_hash = b.eh.specialize(F(sigma))
    g=b.g
    geo = front(model_old,params)
    g.trial(params[-1][0],tight=True)
    b.eh.voronoi_signs(b.eh.arc_data())
    old_norm,_,old_gain,u,w=b.eh.clip_geometry(model_old["parameters"])
    x=b.prev.point_on_faces(old_norm,geo["relaxed"])
    t=geo["t"]
    assert all(g.cross(g.minus(q,p),g.minus(x,p))>0 for p,q in ((u,w),(w,t),(t,u)))
    det=g.cross(*old_norm)
    assert all(g.cross(m,old_norm[1])/det>0 and g.cross(old_norm[0],m)/det>0
               for m in geo["normals"][:2])
    assert g.cross(g.minus(t,u),g.minus(x,u))<0

    weight_data = json.loads(weights.read_bytes())
    keys = [*model_old["ordered_branches"],"H0","0H","1H","HH",*(f"U{i}" for i in range(len(params)))]
    assert weight_data["scale"]=="0.997" and set(keys)==set(weight_data["branches"]) and len(keys)==17
    bounds = {}
    for key in keys:
        p=sign_case(C,model_old,geo,params,key)
        cert=weight_data["branches"][key]
        b.replay_cuts(p,cert["adaptive_cuts"])
        val=p.validate(b.eh.phi_objective(p,b.RIGHT,1),cert["dual"])-arb(3)/2
        assert val<b.arbq(b.F_LIMIT)
        bounds[key]=val
    p=C(params[-1][0])
    b.eh.add_joint_query(p)
    p.linear(b.eh.phi_objective(p,b.RIGHT,1),arb(3)/2+b.arbq(b.F_LIMIT))
    replay=weight_data["pose"]
    b.replay_cuts(p,replay["adaptive_cuts"])
    with redirect_stdout(sys.stderr):
        pose=b.eh.certify(p,b.RIGHT,10,replay,eh_relocation=True,joint_eh=True,revisit_ah=True)["summary"]
    assert arb(pose["minimum_slack_lower_bound"])>0 and arb(pose["F_right_upper"])<0
    new=geometries(sigma,model_old["parameters"],params,g)
    assert new["facets"]==geo["active"] and len(new["facets"])==11
    gain=new["new_eh"]-new["old_eh"]
    assert gain.overlaps(old_gain+geo["gain"])
    return dict(sigma=sigma, right_critical_upper=pack(max(bounds.values())),
                right_by_branch={k:pack(v) for k,v in bounds.items()}, placement=pose,
                active_facets=new["facets"], era_geometry_replays=checks,
                depth=params[-1][0], EH_gain=pack(enclosure(gain),35),
                weights_sha256=sha256(weights.read_bytes()).hexdigest(),
                adapted_front_sha256=front_hash, adapted_sign_sha256=sign_hash, model_sha256=model_hash,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=="__main__":
    assert __debug__
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("comp",type=Path)
    ap.add_argument("sigma")
    ap.add_argument("weights",type=Path)
    a=ap.parse_args()
    print(json.dumps(run(a.comp.resolve(),a.sigma,a.weights),separators=(",",":")))
