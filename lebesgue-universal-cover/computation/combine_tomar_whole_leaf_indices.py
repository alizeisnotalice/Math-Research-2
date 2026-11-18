"""Combine two truncated, certified whole-leaf lower bounds and explicit fixed-point witnesses."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


def sha(data):
    return sha256(data).hexdigest()


def load_shard(directory, row):
    with np.load(directory / f'part{row[0]:04}.npz', allow_pickle=False) as data:
        pos, floor = data['pos'], data['full_floor']
    assert pos.dtype == np.int32 and floor.dtype == np.int64 and pos.shape == floor.shape
    assert np.array_equal(pos, np.unique(pos))
    assert [sha(pos.tobytes()), sha(floor.tobytes())] == row[-2:]
    return pos, floor


def main():
    parser = ArgumentParser(description=__doc__)
    for name in ('library', 'first_report', 'first', 'second_report', 'second',
                 'witness_receipt', 'witness_npz', 'emit'):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    sys.path.insert(0, str(a.library))
    from audit_tomar_multi_target_scores import GAP_EDGES, support
    from audit_tomar_outer_083730 import AREA
    reports = [json.loads(f.read_bytes()) for f in (a.first_report, a.second_report)]
    lhs, rhs = reports
    for key in ('max_index_target', 'canonical_leaves', 'area_denominator', 'coord_denominator', 'gap_edges'):
        assert lhs[key] == rhs[key]
    assert lhs['area_denominator'] == AREA == 2**57
    for key in ('prior_radial', 'manifest', 'root', 'canonical_slice', 'stage_index', 'bucket_index'):
        assert lhs['provenance'][key] == rhs['provenance'][key]
    assert [r[:2] for r in lhs['rows']] == [r[:2] for r in rhs['rows']]
    assert [r[0] for r in lhs['rows']] == sorted(set(r[0] for r in lhs['rows']))
    targets = sorted({F(e['target']) for e in lhs['curve']} |
                     {F('.83742'), F('.83743'), F('.83744')})
    assert max(targets) == F(lhs['max_index_target'])
    goals = [int(t*AREA)+1 for t in targets]
    receipt = json.loads(a.witness_receipt.read_bytes())
    assert receipt['sha256'] == sha(a.witness_npz.read_bytes())
    with np.load(a.witness_npz, allow_pickle=False) as data:
        hints = data['rows']
    assert hints.dtype == np.int64 and hints.shape == (receipt['candidates'], 5)
    assert len(set(map(tuple, hints[:, :2]))) == len(hints)
    a.emit.mkdir(parents=True, exist_ok=True)
    scores, rows = [], []
    subgoal_counts = np.zeros((len(goals), 4), dtype=np.int64)
    for lrow, rrow in zip(lhs['rows'], rhs['rows']):
        part, canonical = lrow[:2]
        x, fx = load_shard(a.first, lrow)
        y, fy = load_shard(a.second, rrow)
        assert np.all(fx < goals[-1]) and np.all(fy < goals[-1])
        z, li, ri = np.intersect1d(x, y, assume_unique=True, return_indices=True)
        common = np.maximum(fx[li], fy[ri])
        bound = common.copy()
        for _, pos, _, center, witness in hints[hints[:, 0] == part]:
            i = np.searchsorted(z, pos)
            if i < len(z) and z[i] == pos:
                bound[i] = max(bound[i], center, witness)
        for j, goal in enumerate(goals):
            subgoal_counts[j] += [(fx < goal).sum(), (fy < goal).sum(),
                                  (common < goal).sum(), (bound < goal).sum()]
        keep = bound < goals[-1]
        z, bound = z[keep], bound[keep]
        np.savez_compressed(a.emit / f'part{part:04}.npz', pos=z, full_floor=bound)
        scores.append(bound)
        rows.append([part, canonical, *[int((bound < g).sum()) for g in goals],
                     sha(z.tobytes()), sha(bound.tobytes())])
    scores = np.sort(np.concatenate(scores))
    np.save(a.emit/'sorted-floor.npy', scores, allow_pickle=False)
    curve = []
    for t, g, counts in zip(targets, goals, subgoal_counts):
        size = int(np.searchsorted(scores, g))
        assert size == counts[3]
        cuts = [size-int(np.searchsorted(scores, g-int(edge*AREA))) for edge in GAP_EDGES]
        worst = int(g-scores[0]) if size else 0
        curve.append(dict(target=str(t), goal=g, unresolved=size, worst_gap_int=worst,
                          radius_surplus_grid_1e8=support(t, worst),
                          gap_counts=np.diff([0, *cuts, size]).tolist(),
                          unresolved_first_second_common_combined=counts.tolist()))
    report = {k:lhs[k] for k in ('area_denominator', 'coord_denominator', 'max_index_target',
                                   'canonical_leaves', 'gap_edges')}
    report.update(sorted_floor_sha256=sha(scores.tobytes()), curve=curve, rows=rows,
                  provenance={**lhs['provenance'], 'first_index':sha(a.first_report.read_bytes()),
                              'second_index':sha(a.second_report.read_bytes()),
                              'integer_witnesses':sha(a.witness_receipt.read_bytes())})
    print(json.dumps(report, separators=(',', ':')))


if __name__ == '__main__':
    main()
