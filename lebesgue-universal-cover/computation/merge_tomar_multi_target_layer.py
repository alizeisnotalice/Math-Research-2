"""Extend a reusable integer index using a completed and separately audited layer run."""
from argparse import ArgumentParser
from fractions import Fraction as F
from io import BytesIO
import json
from pathlib import Path

import numpy as np

from audit_tomar_multi_target_scores import GAP_EDGES, support
from audit_tomar_outer_083730 import AREA, sha


def merge(a):
    base_raw = a.report.read_bytes()
    base = json.loads(base_raw)
    stage_raw = (a.stage/'stage-ledger.json').read_bytes()
    stage = json.loads(stage_raw)
    audit_raw = a.stage_audit.read_bytes()
    audit = json.loads(audit_raw)
    seed_raw = a.seed_report.read_bytes()
    seed = json.loads(seed_raw)
    assert seed['provenance']['whole_leaf_index'] == sha(base_raw)
    assert stage['complete'] and len(stage['audits']) == 4
    assert stage['seed_receipt_sha'] == audit['first_seed_sha256'] == sha(seed_raw)
    assert audit['stage_sha256'] == sha(stage_raw)
    assert F(stage['target']) <= F(base['max_index_target'])
    assert int(F(stage['target'])*AREA)+1 == stage['integer_goal'] == seed['target_integer']
    a.emit.mkdir(parents=True, exist_ok=True)

    goals = [r['goal'] for r in base['curve']]
    stage_maps = [dict(layer['work_files']) for layer in stage['audits']]
    all_scores, rows = [], []
    for part, canonical, *entry in base['rows']:
        name = f'part{part:04}.npz'
        with np.load(a.index/name, allow_pickle=False) as current:
            pos, old = current['pos'], current['full_floor']
        assert [sha(pos.tobytes()), sha(old.tobytes())] == entry[-2:]
        assert pos.dtype == np.int32 and old.dtype == np.int64 and pos.shape == old.shape
        new = old.copy()
        active = pos[old < stage['integer_goal']]
        for j, receipts in enumerate(stage_maps, 1):
            pkt = (a.stage/f'stage{j}'/name).read_bytes()
            assert sha(pkt) == receipts[part]
            with np.load(BytesIO(pkt), allow_pickle=False) as layer:
                passed, pending, lower = layer['accepted'], layer['remaining'], layer['full_floor']
            assert passed.dtype == pending.dtype == np.int32 and lower.dtype == np.int64
            assert lower.shape == passed.shape and np.all(lower >= stage['integer_goal'])
            assert np.array_equal(np.sort(np.r_[passed, pending]), active)
            i = np.searchsorted(pos, passed)
            assert np.all(i < len(pos)) and np.array_equal(pos[i], passed)
            new[i] = lower
            active = pending
        assert active.tolist() == audit['retained'].get(str(part), [])
        keep = new < goals[-1]
        pos, new = pos[keep], new[keep]
        np.savez_compressed(a.emit/name, pos=pos, full_floor=new)
        all_scores.append(new)
        rows.append([part, canonical, *[int((new < g).sum()) for g in goals],
                     sha(pos.tobytes()), sha(new.tobytes())])

    scores = np.sort(np.concatenate(all_scores))
    np.save(a.emit/'sorted-floor.npy', scores, allow_pickle=False)
    curve = []
    for record in base['curve']:
        t, g = F(record['target']), record['goal']
        size = int(np.searchsorted(scores, g))
        cuts = [size-int(np.searchsorted(scores, g-int(edge*AREA))) for edge in GAP_EDGES]
        worst = int(g-scores[0]) if size else 0
        curve.append(dict(target=str(t), goal=g, unresolved=size, worst_gap_int=worst,
                          radius_surplus_grid_1e8=support(t, worst),
                          gap_counts=np.diff([0, *cuts, size]).tolist()))
    return dict(area_denominator=AREA, coord_denominator=base['coord_denominator'],
                max_index_target=base['max_index_target'], canonical_leaves=base['canonical_leaves'],
                sorted_floor_sha256=sha(scores.tobytes()), curve=curve, rows=rows,
                gap_edges=base['gap_edges'],
                provenance={**base['provenance'], 'previous_index':sha(base_raw),
                            'added_stage_ledger':sha(stage_raw), 'added_stage_reuse_audit':sha(audit_raw),
                            'added_stage_seed_report':sha(seed_raw)})


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    for name in ('report', 'index', 'stage', 'seed_report', 'stage_audit', 'emit'):
        parser.add_argument(name, type=Path)
    print(json.dumps(merge(parser.parse_args()), separators=(',', ':')))
