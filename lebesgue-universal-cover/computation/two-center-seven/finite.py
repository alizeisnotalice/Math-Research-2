"""Recompute the finite tangent interpolation for the exact two-centre family."""
from argparse import ArgumentParser
from fractions import Fraction as F
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from zipfile import ZipFile

from flint import arb,ctx

HERE=Path(__file__).resolve().parent
GEOMETRY = HERE.parent / "directional-seven" / "tangent-points.zip"
if sha256(GEOMETRY.read_bytes()).hexdigest() != "e2a48ae7ad4d57d6c66185037024aee895558f879d578b31267c339dcb7f233a":
    raise ValueError("unexpected tangent geometry source")
sys.path.insert(0, str(GEOMETRY))
import certify as t

SC=F(".999999083")
S=F(".99999909")
G=F(120696470013529293,2**57)
H=G/S**2
H_CURVE=G/SC**2


def calibrator_raw(clouds, report, positive=True):
    d=json.loads(report.read_bytes())
    original=clouds["P66" if "P" in report.name else "physicalQ"]
    assert d["denominator"]==original["denominator"]
    assert d["physical_points"]==original["points"]
    assert F(d["eps"])==1-S and F(d["curve_eps"])==1-SC
    assert not positive or (d["diameter_ok"] and d["rows_ok"])


def closest_disk(target):
    lo,hi=3,1000000
    tar=t.ball(target)
    assert t.ball(H_CURVE)*t.disk_factor(lo)**2<tar<t.ball(H_CURVE)*t.disk_factor(hi)**2
    while hi-lo>1:
        n=(lo+hi)//2
        bound=t.ball(H_CURVE)*t.disk_factor(n)**2
        if bound>tar:
            hi=n
        else:
            assert bound<tar
            lo=n
    return [hi,hi+int(hi%2==0)]


def main():
    p=ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    p.add_argument("--add-target",action="append",choices=(".837501",".8375015",".83750152"),default=[])
    a=p.parse_args()
    ctx.prec=256
    a.out.mkdir(parents=True, exist_ok=False)
    inputs = a.out / "source"
    archive = HERE / "source.zip"
    if sha256(archive.read_bytes()).hexdigest() != "17dd03bfeb5177285d9f5afe9e18eecd2778e5f96b1f823703b8506c4b4c6673":
        raise ValueError("unexpected two-centre source")
    with ZipFile(archive) as bundle:
        bundle.extractall(inputs)
    raw=(inputs/"points.json").read_bytes()
    d=json.loads(raw)
    for name in ("nativeP910","strongP910","nativeQ910","expandedQ910","strongQ910"):
        calibrator_raw(d,inputs/(name+".json"))
    targets=(".83750",".837501",".8375015",".83750152")
    tiers={}
    for target in targets:
        k0,k=closest_disk(F(target))
        n0=int((arb.pi()/(2*(t.ball(F(target)/H_CURVE).sqrt()).acos())).ceil().unique_fmpz())
        grids=[k]+[(n0+n-1)//n for n in [3,5,7,9]]
        tiers[target]=dict(minimum_possible_cap=k0, first_odd=k,
                           curve_grids=grids, expected_total=k+sum(n*(grid+1) for n,grid in zip([3,5,7,9],grids[1:]))+77)
    fixed=[]
    for name,c in d.items():
        D=c["denominator"]*S.numerator
        pts=[[x*S.denominator,y*S.denominator] for x,y in c["points"]]
        assert max(t.square_distance(a,b) for a,b in combinations(pts,2)) < D*D
        fixed.append(dict(body=name,points=pts,denominator=D))
    for target in (".83750",*a.add_target):
        r=tiers[target]
        r["curve_certificates"]=[]
        cloud=[]
        for n,k in zip(t.SHAPES,r["curve_grids"],strict=True):
            item, rec=t.build_curve(t.Curve(n),k)
            cloud.append(item)
            r["curve_certificates"].append(rec)
        beta=min(F(1),*(F(c["included_scale"])*S/SC for c in r["curve_certificates"]))
        result=H*beta**2
        assert result>F(target)
        cloud.extend(fixed)
        r.update(actual_counts=[len(g["points"]) for g in cloud], beta=str(beta), area_lower=str(result),
                 decimal_lower=t.decimal_floor(result), max_count=max(len(g["points"]) for g in cloud))
        r["total_count"]=sum(r["actual_counts"])
        assert r["expected_total"]==r["total_count"]
        canonical=t.packed(cloud)
        r["point_coordinates_sha256"]=sha256(canonical).hexdigest()
        name="tangent-910-points.json.gz" if target==".83750" else f"tangent-910-points-{target[1:]}.json.gz"
        (a.out/name).write_bytes(gzip.compress(canonical,compresslevel=8,mtime=0))
        print(target,r["actual_counts"],r["decimal_lower"],flush=True)
    out=dict(model_statement="five curves at s_c and two rational point groups with diameter at most s",
             scale_curves=str(SC),final_normalizer=str(S),area_unnormalized=str(G),area_normalized=str(H),
             source_points_sha256=sha256(raw).hexdigest(), tiers=tiers, precision=ctx.prec)
    (a.out/"finite-thresholds.json").write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    summary={k:[v["minimum_possible_cap"],v["first_odd"],v["expected_total"]] for k,v in tiers.items()}
    print(summary,flush=True)


if __name__=="__main__":
    main()
