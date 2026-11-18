"""Certify two fixed X18 poses throughout the path and its Hausdorff tube."""
import json
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx
from certify_six_anchor_path import FIVE_SPEED, ar, data_files, midpoint, motion, pose_at
from certify_x18_lower_hosts import inner_hull, positivity


def verify(directory):
    ctx.prec = 192
    inputs = data_files(directory)
    start = list(map(F, inputs["free-p30-reoptimized-pose.json"]))
    end = list(map(F, inputs["structural-six-escape-input.json"]["pose"]))
    motion(start, end, inputs["finite-point-leaf-input.json"], arb.pi())
    poses = json.loads((directory / "x18-anchor-path-poses.json").read_bytes())
    cloud = json.loads((directory / "widened-projected-obstruction-18-points.json").read_bytes())
    assert [pose["sign"] for pose in poses] == [1, -1] and len(cloud["points"]) == 18
    assert max((a-c)**2+(b-d)**2 for a,b in cloud["points"] for c,d in cloud["points"]) < cloud["denominator"]**2
    cells = 16
    gaps = {1: [], -1: []}
    for index in range(cells):
        center = pose_at(start[:11], end[:11], midpoint(index, cells))
        polygon, denominator = inner_hull(list(map(ar, center)), 128)
        for pose in poses:
            gaps[pose["sign"]].append(positivity(polygon, denominator, pose, cloud))
    margins = {sign: min(row) - FIVE_SPEED/(2*cells) for sign, row in gaps.items()}
    assert min(margins.values()) > 0
    return {"poses":poses, "midpoint_margins":{k:list(map(str,v)) for k,v in gaps.items()},
            "whole_path_clearance":{k:str(v) for k,v in margins.items()},
            "display":{k:float(v) for k,v in margins.items()}}


if __name__ == "__main__":
    print(json.dumps(verify(Path(__file__).resolve().parent), indent=2))
