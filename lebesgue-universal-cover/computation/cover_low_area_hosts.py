"""Cover convex hosts of five fixed test bodies by explicitly placed center hulls."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

import numpy as np
from flint import arb, ctx
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST

HERE = Path(__file__).resolve().parent
GRID = 2**46
TARGET = F(".83732")
FIVE = F(120666205824033364, 2**57)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def center_boxes(tree):
    manifest = tree / "SHA256SUMS"
    assert digest(manifest) == PINNED_MANIFEST
    pins = {name.removeprefix("./"):h for h,name in (line.split() for line in manifest.read_text().splitlines())}
    chart = HERE / "tomar-083732-retained.json"
    assert digest(chart) == "605e3b3b2d6268a41faba14300c61dcebfd27c726d549e13916c08261bc855a9"
    ids, centers, widths = [], [], []
    for part, positions in json.loads(chart.read_bytes()).items():
        header = tree / f"part{int(part):04}.json"
        codefile = header.with_suffix(".npy")
        assert all(digest(p) == pins[p.name] for p in [header, codefile])
        meta = json.loads(header.read_bytes())
        code = np.load(codefile, allow_pickle=False)
        prefix = np.r_[1, 1 + np.cumsum(np.where(code < 11, 1, -1))]
        for pos in positions:
            point, width = leaf_box(prefix, code, int(part), pos, meta)
            for values, dest in [(point, centers), (width, widths)]:
                row = [F(v)*GRID for v in values]
                assert all(x.denominator == 1 for x in row)
                dest.append(list(map(int,row)))
            ids.append([int(part),pos])
    assert len(ids) == 705
    return ids, np.asarray(centers), np.asarray(widths)


def atlas(tree):
    ids, centers, widths = center_boxes(tree)
    ctx.prec = 160
    # For a circular arc the rotational derivative of its support is the
    # tangential component of the center. At a vertex it is that of the vertex.
    # In both normal sectors the absolute angular offset is at most pi/(2n).
    for n in (5,7,9):
        assert (arb.pi()/(2*n)).tan()/2 < arb(163)/1000
    radius = arb(703)/1000
    spike = (arb.pi()-(1/(2*radius)).acos())/4 + (radius**2-arb(1)/4).sqrt()/2
    assert spike > arb(TARGET.numerator)/TARGET.denominator
    gap = TARGET-FIVE
    delta = F(548,100000)
    assert gap**2*(4*F(".703")**2-1) < delta**4

    size = len(ids)
    costs = np.zeros((size,size))
    placements = centers / GRID
    for i, j, angle in [(0,1,None),(3,4,2),(6,7,5),(9,10,8)]:
        def offset(axis):
            return abs(placements[:,None,axis]-placements[None,:,axis]) + widths[None,:,axis]/GRID
        block = np.hypot(offset(i),offset(j))
        if angle is not None:
            block += .163*offset(angle)
        costs = np.maximum(costs,block)
    # Doubles only select the representatives and assignments.
    def verify(rep, leaf, threshold):
        changes = abs(centers[rep]-centers[leaf])+widths[leaf]
        for i,j,angle in [(0,1,None),(3,4,2),(6,7,5),(9,10,8)]:
            x,y = map(int, (changes[i],changes[j]))
            angular = F(163*int(changes[angle]),1000) if angle is not None else 0
            budget = threshold*GRID - angular
            assert budget >= 0 and x*x+y*y <= budget**2

    levels = []
    for threshold in map(F, (".10",".06",".04",".03",".02",".015",".012",".01")):
        cover = costs < float(threshold)
        if not cover.any(axis=0).all():
            continue
        live = np.ones(size,dtype=bool)
        groups = []
        while live.any():
            rep = int(np.argmax(cover[:,live].sum(axis=1)))
            chosen = np.flatnonzero(live & cover[rep])
            assert len(chosen)
            for leaf in chosen:
                verify(rep, int(leaf), threshold)
            live[chosen] = False
            groups.append(dict(center_numerators=centers[rep].tolist(), representative=ids[rep],
                               leaves=[ids[i] for i in chosen]))
        assert sum(len(g["leaves"]) for g in groups)==size
        levels.append(dict(host_radius=str(threshold+delta), hull_radius=str(threshold),
                           count=len(groups), groups=groups))
    return dict(target=str(TARGET), five_area_floor=str(FIVE), host_to_own_hull=str(delta),
                coordinate_denominator=GRID, levels=levels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree",type=Path)
    args = parser.parse_args()
    print(json.dumps(atlas(args.tree),separators=(",",":")))
