"""Recompute exact four-dimensional union by Cartesian box intersections."""
import hashlib
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction
from functools import cache
from pathlib import Path

computation = Path(__file__).resolve().parent
prior_raw = (computation / '027-b3-pruning.json').read_bytes()
current_raw = (computation / '029-four-dimensional-pruning.json').read_bytes()
prior = json.loads(prior_raw)
current = json.loads(current_raw)
assert prior["tau"] == current["tau"] == "419/500"
assert prior["coordinate_denominator"] == current["coordinate_denominator"]
assert prior["source_global_box"]["sha256"] == current["source_global_box"]["sha256"]
prior_hash = hashlib.sha256(prior_raw).hexdigest()
assert current["comparison_b3_certificate"]["sha256"] == prior_hash
scale = int(prior["coordinate_denominator"])
b_root = tuple(map(int, prior["translation_domain_scaled"]))
f_root = (-scale * 6 // 25, scale * 6 // 25)
root = (b_root, b_root, f_root, f_root)

def volume(box):
    return math.prod(hi - lo for lo, hi in box)

def intersection_volume(first, second):
    return math.prod(max(0, min(hi, end) - max(lo, start))
                     for (lo, hi), (start, end) in zip(first, second, strict=True))

def parse_prior(leaf):
    return tuple(tuple(map(int, leaf[axis])) for axis in ('x_scaled', 'y_scaled'))

prior_pruned = [parse_prior(leaf) for leaf in prior["leaves"] if leaf["status"] == "prune"]
assert all(intersection_volume(left, right) == 0
           for i, left in enumerate(prior_pruned)
           for right in prior_pruned[i+1:])
assert all(volume(box) == intersection_volume(box, root[:2]) for box in prior_pruned)

@cache
def prior_overlap(projection):
    return sum(intersection_volume(projection, prior_box) for prior_box in prior_pruned)

prior_volume = sum(map(volume, prior_pruned)) * volume(root[2:])
current_pruned = overlap = overlap_from_retained = remainder = 0
residual_by_reason = defaultdict(int)
overlap_by_relation = defaultdict(int)
for leaf in current["leaves"]:
    box = tuple(tuple(map(int, pair)) for pair in leaf["box_scaled"])
    box_volume = volume(box)
    assert box_volume == intersection_volume(box, root)
    in_prior = prior_overlap(box[:2]) * volume(box[2:])
    if leaf["status"] == "prune":
        current_pruned += box_volume
        overlap += in_prior
        overlap_by_relation[leaf["b3_projection_relation"]] += in_prior
    else:
        assert leaf["status"] == "retain"
        overlap_from_retained += in_prior
        remainder += box_volume - in_prior
        residual_by_reason[leaf["reason"]] += box_volume - in_prior
assert overlap + overlap_from_retained == prior_volume
union_volume = current_pruned + overlap_from_retained
assert union_volume + remainder == volume(root)
assert union_volume == prior_volume + current_pruned - overlap

fractions = {
    'only_027': prior_volume - overlap,
    'only_029': current_pruned - overlap,
    'both': overlap,
    '027_total': prior_volume,
    '029_total': current_pruned,
    'union': union_volume,
    'complement': remainder,
}

def fraction_text(values):
    return {label: str(Fraction(v, volume(root))) for label, v in values.items()}

print(json.dumps({
    'python': sys.version,
    '027_sha256': prior_hash,
    '029_sha256': hashlib.sha256(current_raw).hexdigest(),
    'root_volume': str(Fraction(volume(root), scale ** 4)),
    'fraction': fraction_text(fractions),
    'absolute_volume': {label: str(Fraction(v, scale ** 4)) for label, v in fractions.items()},
    'intersection_by_029_relation': fraction_text(overlap_by_relation),
    'complement_by_029_reason': fraction_text(residual_by_reason),
    'unique_projected_rectangles': prior_overlap.cache_info().currsize,
}, indent=2, sort_keys=True))
