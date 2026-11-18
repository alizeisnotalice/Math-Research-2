"""Certify the .83742 shallow-cap transfer on the fixed 13,368 rational points."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("computation","tree","table32","table64","work","cache"):
        parser.add_argument(name,type=Path)
    parser.add_argument("--margin-kit",type=Path)
    args=parser.parse_args()
    assert __debug__
    comp=args.computation.resolve()
    sys.path.insert(0,str(comp))
    from certify_polygon_083731 import ORDERS, ar, arb, ctx, lower_radius
    from certify_tomar_six_high import SourceTree
    from select_tomar_083731_scores import radius_rows
    ctx.prec=160
    work=args.work.resolve()
    work.mkdir(exist_ok=False)
    old_dir=work/"geometry"
    subprocess.run([sys.executable,"-I","-B",comp/"certify_signed_shrink_083741.py",
                    comp,args.tree,args.table32,args.table64,old_dir,args.cache,
                    "--base-kit",comp/"five-affine-selected-083741.zip",
                    "--margin-kit",comp/"five-affine-shallow-cap-34-083742.zip"],
                   check=True,stdout=subprocess.DEVNULL)
    old=read(old_dir/"replay.json")
    assert old==read(comp/"signed-shallow-cap-083741-replay.json") and old["complete"]
    terminal=read(comp/"finite-point-lift-083742.json")
    seed=read(comp/"tomar-083742-seed.json")
    stage=read(comp/"tomar-083742-stage-ledger.json")
    screen=read(comp/"tomar-083742-signed-screen-audit.json")
    scores=read(comp/"tomar-083742-score-layers-audit.json")
    retained=read(comp/"tomar-083742-retained.json")
    T=F(".83742")
    assert all(F(v["target"])==T for v in (terminal,seed,stage,screen,scores))
    assert stage["complete"] and stage["parts"]==1024 and screen["complete"]
    assert screen["geometry_source_sha256"]==stage["source_sha"]==old["signed_kernel_sha256"]
    assert digest(comp/"tomar-083742-seed.json")==stage["seed_receipt_sha"]==scores["first_seed_sha256"]
    assert scores["stage_sha256"]==terminal["outer_audit_sha256"]==digest(comp/"tomar-083742-stage-ledger.json")
    assert scores["retained"]==retained and terminal["retained_sha256"]==digest(comp/"tomar-083742-retained.json")
    assert seed["tree_manifest_sha256"]==digest(args.tree/"SHA256SUMS")
    for key,name in (("weighted_sha256","tomar_six_083730_guard.cpp"),("point_sha256","finite-point-leaf-input.json"),
                     ("bank_sha256","finite-point-global-input.json"),("base_sha256","tomar-six-high-input.json")):
        assert terminal[key]==read(comp/"finite-point-lift-083741.json")[key]==digest(comp/name)
    goal=F(seed["target_integer"],2**57)
    assert int(goal*2**57)==stage["integer_goal"]==screen["goal"]==120684940798483225
    radial=read(comp/"five-body-radial-affine.json")
    assert digest(comp/"five-body-radial-affine.json")==seed["provenance"]["radial_certificate"]
    assert radius_rows(radial,F(seed["radial_target"]))==seed["radii"]
    assert all(F(row["integer_radius"],2**44)<F(".069") and F(row["integer_area_lower"],2**57)>goal for row in seed["radii"])

    previous=read(comp/"shallow-cap-audit.json")
    s,r,R=F(old["curve_scale"]),F(1,2),F(".703")
    assert s==F(".999999035")
    root=SourceTree(args.tree).root
    for j,n in enumerate(ORDERS[1:]):
        k=3*j
        dist=min(F(root["hi"][k]),F(root["hi"][k+1]),-F(root["lo"][k]),-F(root["lo"][k+1]))
        ins=ar(s)*lower_radius(n,arb.pi())
        assert (arb.pi()*(ar(s*r)**2+ins**2)/2+ar(dist)*(ar(s*r)+ins))>ar(goal)
    disk_area=(arb.pi()-(ar(s*r)/ar(R)).acos())*ar(s*r)**2+ar(s*r)*(ar(R)**2-ar(s*r)**2).sqrt()
    assert disk_area>ar(goal) and R**2<2*(s*r)**2
    h=F(".01451")
    max_target=s*s*F(".837295")+s*r*h*h/R
    assert max_target>F(".83744312")>goal
    m=F(previous["rho"])*F(".528619")/(1-s)
    a=r-F(".069")
    v=(1+s)*r*r+m*m*(1-s)
    slope=1+s-(2*m*r+v)/(4*r*r)
    margin=4*(a-m)*s*s-(1+s)*h-max(0,(v-2*m*(r+h))/(2*r+h))
    assert m<a<r and slope>0 and margin>F(".00000895")
    loss=F(old["loss"])
    assert F(".837423")-loss>goal
    pure={tuple(x) for x in terminal["five_only"]}
    ordinary={tuple(x[:2]) for x in terminal["rows"]}
    keys={(int(i),j) for i,items in retained.items() for j in items}
    assert len(ordinary)==2300 and len(pure)==66 and not pure&ordinary and pure|ordinary==keys
    assert len(keys)==stage["audits"][-1]["remaining"]==2366
    factor=F(old["rational_scale"])
    result=dict(complete=False,terminal_sha256=digest(comp/"finite-point-lift-083742.json"),
                old_geometry_sha256=digest(old_dir/"replay.json"),coordinate_sha256=old["coordinate_sha256"],
                signed_kernel_sha256=old["signed_kernel_sha256"],integer_goal=str(goal),cap_height=str(h),
                shallow_comparison=str(margin),comparison_derivative=str(slope),cap_target_limit=str(max_target),
                exceptional_loss=str(loss),rational_scale=str(factor),conditional_floor=str(factor**2*goal))
    if args.margin_kit:
        base=comp/"five-affine-selected-083742-final.zip"
        assert digest(base)==terminal["five_kit_sha256"]
        strong=[]
        with ZipFile(base) as orig, ZipFile(args.margin_kit) as high:
            for i,j in sorted(pure):
                stem=f"expected/part{i}-pos{j}"
                a=json.loads(orig.read(stem+".json"));b=json.loads(high.read(stem+".json"))
                for field in ("root","upstream","original_tree_sha256","original_header_sha256"):
                    assert a[field]==b[field]
                assert (b["part"],b["pos"])==(i,j) and F(b["threshold"])>=F(".837423") and not b["result"]["remaining"]
                assert b["trace_sha256"]==sha256(high.read(stem+".bin")).hexdigest()
                strong.append([i,j,b["threshold"],b["result"]["vis"],b["trace_sha256"]])
        floor=goal*factor**2
        assert floor>F(".83742123227")
        result.update(complete=True,high_margin_sha256=digest(args.margin_kit),high_rows=strong,
                      area_floor=str(floor),decimal_floor=".83742123227")
    (work/"replay.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="high_rows"},indent=2))


if __name__=="__main__":
    main()
