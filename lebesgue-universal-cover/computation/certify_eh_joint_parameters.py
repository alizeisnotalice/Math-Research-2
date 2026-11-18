#!/usr/bin/env python3
"""Specialize the certified EH relocation model to a rational Gibbs angle."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import audit_optimized_eh_chord as geometry
import certify_eh_relocation_adaptive as adaptive
import certify_eh_relocation_root_condition as root_condition
from certify_eh_general_relocation import exact_fraction


def specialize(sigma):
    # All geometric functions reference this one exact angle at call time.
    geometry.SIGMA=exact_fraction(sigma)*geometry.PI/180
    model_path=Path(root_condition.__file__)
    source=model_path.read_text()
    old='sigma=arb("1.5494009336023418793")*PI/180'
    replacement='sigma=geometry.SIGMA'
    assert source.count(old)==1
    specialized=source.replace(old,replacement)
    namespace={"__name__":"specialized_eh_model","__file__":str(model_path),"geometry":geometry}
    exec(compile(specialized,str(model_path),"exec"),namespace)
    return namespace["RootConditionProgram"],hashlib.sha256(specialized.encode()).hexdigest()


def gain_endpoint(ball):
    mantissa,exponent=(int(c) for c in ball.man_exp())
    return Fraction(mantissa)*Fraction(2)**exponent


def run(args):
    sigma=Fraction(args.sigma)
    assert 0<sigma<30
    model_type,model_hash=specialize(sigma)
    geometry_checks=geometry.trial(args.delta,tight=True)
    program=model_type(args.delta)
    adaptive_defaults=adaptive.multipliers

    def multipliers(name):
        if name=="K-C3":
            return tuple(map(Fraction,("0","1.35","1.5","1.65","1.75","2","2.25","2.5","3","5")))
        return adaptive_defaults(name)

    adaptive.multipliers=multipliers
    saved=None
    if args.replay:
        saved=json.loads(args.replay.read_text())
        assert Fraction(saved["summary"]["sigma_degrees"])==sigma
        assert saved["summary"]["delta"]==args.delta
        assert saved["summary"]["specialized_model_sha256"]==model_hash
        assert Fraction(saved["summary"]["angle_root_range_degrees"][1])==Fraction(args.right)
        for index,direction in saved["adaptive_cuts"]:
            program.cut(index,direction,record=True)
    result=adaptive.certify(program,Fraction(args.right),args.max_depth,saved)
    summary=result["summary"]
    summary["sigma_degrees"]=str(sigma)
    summary["specialized_model_sha256"]=model_hash
    summary["certificate_wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    gain=geometry.arb(geometry_checks["gain_area"])
    summary["cut_gain_area_interval_exact"]=[str(gain_endpoint(v)) for v in (gain.lower(),gain.upper())]
    if args.old_area:
        old=json.loads(args.old_area.read_text())
        assert Fraction(old["sigma_degrees"])==sigma
        a,b=map(Fraction,old["area_interval"])
        low,high=map(Fraction,summary["cut_gain_area_interval_exact"])
        summary["candidate_area_interval_exact"]=[str(a-high),str(b-low)]
        summary["old_area_input_sha256"]=hashlib.sha256(args.old_area.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sigma")
    parser.add_argument("delta")
    parser.add_argument("right")
    parser.add_argument("output",type=Path)
    parser.add_argument("--max-depth",type=int,default=8)
    parser.add_argument("--old-area",type=Path)
    parser.add_argument("--replay",type=Path)
    run(parser.parse_args())
