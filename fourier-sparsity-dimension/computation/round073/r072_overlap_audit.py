"""All 374 five-bit kernels, with exact active and minimum selected-coset overlaps."""
from collections import Counter
import json
import platform
from r072_all_quotients import affine_hull, linear_subspaces


def walsh_graph():
    # Same two truth tables a=1, b=8 as in the main audit, omitting a shared factor.
    out = {}
    for u in range(4):
        for v in range(4):
            a, b = int(u == 0), int(v == 3)
            freq = (v << 3) | (u << 1) | (a ^ b)
            out[freq] = 1 - 2 * (((u & v).bit_count() + a) % 2)
    return out


def masks_to_positions(mask):
    return tuple(i for i in range(mask.bit_length()) if mask & (1 << i))


def feasible(options, size, cap):
    loads = [0]*size
    def choose(i):
        if i == len(options):
            return True
        for o in options[i]:
            if any(loads[x] >= cap for x in o):
                continue
            for x in o:
                loads[x] += 1
            found = choose(i+1)
            for x in o:
                loads[x] -= 1
            if found:
                return True
        return False
    return choose(0)


def overlaps(H, signs):
    restrictions = {}
    for x in range(32):
        key = tuple((x & h).bit_count() % 2 for h in sorted(H))
        restrictions.setdefault(key, x)
    physical = list(restrictions.values())
    bands = {}
    for a in signs:
        bands.setdefault(min(a ^ z for z in H), []).append(a)
    options = []
    actives = [0]*len(physical)
    filling = 0
    for points in bands.values():
        anchor = points[0]
        diffs = [z ^ anchor for z in points]
        hull_size = len(affine_hull(points))
        filling += hull_size
        masks = {}
        for i, x in enumerate(physical):
            phase = tuple((x & z).bit_count() % 2 for z in diffs)
            value = sum(signs[z]*(-1)**s for z, s in zip(points, phase))
            if value:
                actives[i] += 1
                masks[phase] = masks.get(phase, 0) | (1 << i)
        assert masks
        assert all(m.bit_count()*hull_size == len(physical) for m in masks.values())
        options.append([masks_to_positions(m) for m in masks.values()])
    options.sort(key=lambda op: (len(op), -len(op[0])))
    lower = (sum(len(op[0]) for op in options)+len(physical)-1)//len(physical)
    actual = max(actives)
    optimal = next(x for x in range(lower, actual+1) if feasible(options, len(physical), x))
    assert optimal*filling >= len(bands)**2
    return optimal, actual, filling


signs = walsh_graph()
hist = Counter()
strict = 0
for H in linear_subspaces(5):
    opt, true, filling = overlaps(H, signs)
    assert opt <= true
    strict += opt < true
    hist[(opt, true)] += 1
print(json.dumps({"python": platform.python_version(), "seed_truth_table_masks": [1, 8],
                  "kernels": sum(hist.values()), "strict_overlap_improvements": strict,
                  "overlap_histogram": [{"optimal": a, "true": b, "count": count}
                                        for (a, b), count in sorted(hist.items())]}, indent=2))
