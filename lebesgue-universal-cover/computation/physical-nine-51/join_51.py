"""Combine the proved .83751 exterior and disjoint exact original core forests."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
for name in ("old_core", "exterior", "p128", "old_q", "new_q", "p15", "output"):
    p.add_argument(name, type=Path)
p.add_argument("core_forests", nargs="+", type=Path)
a = p.parse_args()
sha = lambda q: sha256(q.read_bytes()).hexdigest()
read = lambda q: json.loads(q.read_bytes())
target = Fraction(".83751")
area = 1 << 57
grid_floor = int(target * area) + 1
roots = {f"{part}/{pos}" for part, positions in read(a.old_core).items() for pos in positions}
exterior = read(a.exterior)
assert len(roots) == exterior["open_originals_after_this_proof"] == 5442
assert sha(a.old_core) == exterior["historical_retained_sha256"]
assert Fraction(exterior["certified_minimum_outside_5442"]) == Fraction(grid_floor, area)
assert exterior["canonical_original_leaves"] == exterior["excluded_canonical_exterior"] + len(roots)

families = ("P128", "oldQ", "newQ", "P15")
points = {}
for family, path, size in zip(families, (a.p128,a.old_q,a.new_q,a.p15), (128,11,11,15), strict=True):
    doc = read(path)
    coords = doc["points"]
    denominator = doc.get("denominator", doc.get("point_denominator"))
    assert len(coords) == size and isinstance(denominator,int) and denominator > 0
    assert all(type(q) is int for pair in coords for q in pair)
    numerator = max((x-u)**2+(y-v)**2 for (x,y),(u,v) in combinations(coords,2))
    assert numerator < denominator**2
    points[family] = dict(point_count=size, max_squared_diameter=str(Fraction(numerator,denominator**2)), sha256=sha(path))

seen, pending, families_used, sources = set(),set(),Counter(),[]
n_closed = 0
nodes = 0
known = {}
for file in a.core_forests:
    data = read(file)
    digest = sha(file)
    assert Fraction(data["target"]) == target
    closed = data["strict_selected_roots"]
    waiting = set(data["pending"])
    names = set(closed) | waiting
    assert not (set(closed) & waiting) and len(names) == data["original_roots"]
    assert len(closed) == data["closed"] and names <= roots and not (seen & names)
    seen.update(names)
    pending.update(waiting)
    n_closed += len(closed)
    nodes += data["all_geometric_visits"]
    for k,r in closed.items():
        assert r["calls"] == 2*(r["five"]+sum(r[label] for label in families))-1
        need = [label for label in families if r[label]]
        families_used["+".join(need) if need else "pure-five"] += 1
        known[k] = dict(families=need, visits=r["calls"], source_sha256=digest)
    sources.append(dict(original_roots=len(names), closed=len(closed), sha256=digest))
assert seen == roots
assert n_closed + len(pending) == len(roots)
common = dict(strict_target=str(target), integer_area_numerator=grid_floor, integer_area_denominator=area,
              complete_global=not pending, old_core_closed=n_closed, old_core_pending=sorted(pending),
              canonical_leaves=exterior["canonical_original_leaves"], exterior_closed=exterior["excluded_canonical_exterior"],
              closed_core_families=dict(sorted(families_used.items())), finite_inputs=points,
              total_core_geometric_visits=nodes, exterior_audit_sha256=sha(a.exterior), core_sources=sources,
              known_originals={k:r for k,r in sorted(known.items())})
if not pending:
    common["universal_cover_lower_bound"] = str(Fraction(grid_floor,area))
a.output.write_text(json.dumps(common,indent=2)+"\n")
print(json.dumps({k: common[k] for k in ("complete_global","old_core_closed","total_core_geometric_visits")}))
