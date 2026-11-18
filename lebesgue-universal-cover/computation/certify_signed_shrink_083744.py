"""Apply the checked shallow-cap comparison to the .83744 signed cover, dense point leaves and five-only exceptions."""
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
    parser.add_argument("--special-source",type=Path)
    parser.add_argument("--terminal",type=Path)
    parser.add_argument("--five-kit",type=Path)
    parser.add_argument("--margin-kit",type=Path)
    args=parser.parse_args()
    assert __debug__ and (args.terminal is None)==(args.five_kit is None)
    assert args.margin_kit is None or args.terminal is not None
    comp=args.computation.resolve()
    sys.path.insert(0,str(comp))
    from certify_polygon_083731 import ORDERS, ar, arb, ctx, lower_radius
    from certify_tomar_six_high import SourceTree
    from select_tomar_083731_scores import radius_rows
    ctx.prec=160
    work=args.work.resolve();work.mkdir(exist_ok=False)
    olddir=work/"previous"
    subprocess.run([sys.executable,"-I","-B",comp/"certify_signed_shrink_083743.py",comp,args.tree,
                    args.table32,args.table64,olddir,args.cache,"--terminal",comp/"finite-point-lift-083743.json","--five-kit",comp/"five-affine-selected-083743-final.zip"],
                   stdout=subprocess.DEVNULL,check=True)
    old=read(olddir/"replay.json")
    assert old["complete"] and old==read(comp/"signed-shallow-cap-083743-replay.json")
    front="tomar-083744-"
    seed=read(comp/(front+"seed.json"));stage=read(comp/(front+"stage-ledger.json"))
    screen=read(comp/(front+"signed-screen-audit.json"));scores=read(comp/(front+"score-layers-audit.json"))
    roots=read(comp/(front+"retained.json"))
    assert all(F(z["target"])==F(".83744") for z in (seed,stage,screen,scores))
    assert stage["complete"] and stage["parts"]==1024 and screen["complete"]
    assert screen["geometry_source_sha256"]==stage["source_sha"]==read(comp/"signed-shallow-cap-083742-replay.json")["signed_kernel_sha256"]
    assert digest(comp/(front+"seed.json"))==stage["seed_receipt_sha"]==scores["first_seed_sha256"]
    assert scores["stage_sha256"]==digest(comp/(front+"stage-ledger.json")) and roots==scores["retained"]
    assert seed["tree_manifest_sha256"]==digest(args.tree/"SHA256SUMS")
    goal=F(stage["integer_goal"],2**57)
    assert int(goal*2**57)==seed["target_integer"]==screen["goal"]==120687823102244742
    radial=read(comp/"five-body-radial-affine.json")
    assert digest(comp/"five-body-radial-affine.json")==seed["provenance"]["radial_certificate"]
    assert radius_rows(radial,F(seed["radial_target"]))==seed["radii"]
    assert all(F(row["integer_radius"],2**44)<F(".069") and F(row["integer_area_lower"],2**57)>goal for row in seed["radii"])
    s,r,R=F(".999999035"),F(1,2),F(".703")
    root=SourceTree(args.tree).root
    for i,n in enumerate(ORDERS[1:]):
        k=3*i
        d=min(F(root["hi"][k]),F(root["hi"][k+1]),-F(root["lo"][k]),-F(root["lo"][k+1]))
        ins=ar(s)*lower_radius(n,arb.pi())
        assert arb.pi()*(ar(s*r)**2+ins**2)/2+ar(d)*(ar(s*r)+ins)>ar(goal)
    assert (arb.pi()-(ar(s*r)/ar(R)).acos())*ar(s*r)**2+ar(s*r)*(ar(R)**2-ar(s*r)**2).sqrt()>ar(goal)
    assert F(old["allowed_cap_target"])>F(".83744312")>goal
    loss=F(old["steiner_loss"]);high=F(".837445")
    assert high-loss>goal
    factor=F(read(comp/"signed-shallow-cap-083742-replay.json")["rational_scale"]);floor=factor*factor*goal
    assert floor>F(".83744123230")
    keys={(int(k),i) for k,ids in roots.items() for i in ids}
    assert len(keys)==stage["audits"][-1]["remaining"]==2963
    out=dict(complete=False,goal=str(goal),conditional_floor=str(floor),coordinate_sha256=old["coordinate_sha256"],
             unclipped_source_receipt=digest(olddir/"replay.json"),allowed_cap_target=old["allowed_cap_target"],
             steiner_loss=str(loss),margin_at_0837445=str(high-loss-goal),retained_count=len(keys))
    if args.terminal:
        terminal=read(args.terminal);prior=read(comp/"finite-point-lift-083743.json")
        assert F(terminal["target"])==F(".83744") and terminal["retained_sha256"]==digest(comp/(front+"retained.json"))
        assert terminal["outer_audit_sha256"]==scores["stage_sha256"]
        for k in ("base_sha256","point_sha256","bank_sha256","weighted_sha256"):
            assert terminal[k]==prior[k]
        pure=set(map(tuple,terminal["five_only"]));regular={tuple(z[:2]) for z in terminal["rows"]}
        special_comp=(args.special_source or comp).resolve()
        dense=terminal["dense"]
        dense_path=special_comp/dense["filename"]
        assert digest(dense_path)==dense["sha256"]
        dense_summary=json.loads(subprocess.check_output([sys.executable,"-I","-B",
                special_comp/"check_scaled_p30_leaves_083744.py",comp,args.tree,dense_path]))
        assert dense_summary["complete"] and F(dense_summary["material_scale"])<=s
        dense_keys=set(map(tuple,read(dense_path)["original_leaves"]))
        assert dense_keys==set(map(tuple,dense["original_leaves"])) and len(dense_keys)==dense_summary["leaves"]
        assert (180,24225) in dense_keys
        groups=(pure,regular,dense_keys)
        assert set.union(*groups)==keys and sum(map(len,groups))==len(keys) and len(regular)==len(terminal["rows"])
        out["mesh_leaves"]=len(dense_keys);out["scaled_mesh_sha256"]=dense["sha256"]
        assert digest(args.five_kit)==terminal["five_kit_sha256"]
        strength=args.margin_kit or args.five_kit
        checked=[]
        with ZipFile(args.five_kit) as base,ZipFile(strength) as kit:
            for i,j in sorted(pure):
                stem=f"expected/part{i}-pos{j}"
                a=json.loads(base.read(stem+".json"));b=json.loads(kit.read(stem+".json"))
                assert all(a[k]==b[k] for k in ("root","original_header_sha256","original_tree_sha256","upstream"))
                assert (i,j)==(b["part"],b["pos"]) and not b["result"]["remaining"]
                assert b["trace_sha256"]==sha256(kit.read(stem+".bin")).hexdigest()
                assert F(b["threshold"])-loss>goal
                checked.append([i,j,b["threshold"],b["result"]["vis"],b["trace_sha256"]])
        out.update(complete=True,terminal_sha256=digest(args.terminal),high_pack_sha256=digest(strength),
                   ordinary=len(regular),pure=len(pure),high_rows=checked,area_floor=str(floor),decimal_floor=".83744123230")
    (work/"replay.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="high_rows"},indent=2))


if __name__=="__main__":
    main()
