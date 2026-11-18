"""Turn an authenticated multi-target residue into seeds for an exact four-layer replay."""
from argparse import ArgumentParser
from fractions import Fraction as F
import json
from pathlib import Path

from flint import ctx
import numpy as np

from audit_tomar_outer_083730 import AREA, sha
from select_tomar_083731_scores import radius_rows


def prepare(a):
    ctx.prec = 160
    raw = a.index_report.read_bytes()
    report = json.loads(raw)
    assert report['area_denominator'] == AREA
    assert 0 < a.target <= F(report['max_index_target']) and a.radial_target >= a.target
    radial = a.radial_certificate.read_bytes()
    radial_hash = sha(radial)
    assert radial_hash == report['provenance']['prior_radial']
    radii = radius_rows(json.loads(radial), a.radial_target)
    goal = int(a.target*AREA)+1
    a.emit.mkdir(parents=True, exist_ok=True)
    records = []
    for part, canonical, *row in report['rows']:
        name = f'part{part:04}.npz'
        with np.load(a.index/name, allow_pickle=False) as z:
            pos, lower = z['pos'], z['full_floor']
        assert pos.dtype == np.int32 and lower.dtype == np.int64 and pos.shape == lower.shape
        assert np.array_equal(pos, np.unique(pos))
        assert [sha(pos.tobytes()), sha(lower.tobytes())] == row[-2:]
        keep = lower < goal
        todo, floors = pos[keep], lower[keep]
        np.savez_compressed(a.emit/name, pos=todo, score=floors)
        records.append([part, canonical, canonical-len(todo), 0, len(todo),
                        sha(todo.tobytes()), sha(floors.tobytes())])
    assert sum(r[1] for r in records) == report['canonical_leaves']
    remaining = sum(r[4] for r in records)
    for row in report['curve']:
        if row['goal'] == goal:
            expected = row['unresolved'] if 'unresolved' in row else row['counts'][4]
            assert remaining == expected
    return dict(target=str(a.target), radial_target=str(a.radial_target), target_integer=goal,
                radii=radii, tree_manifest_sha256=report['provenance']['manifest'],
                provenance={'radial_certificate':radial_hash, 'whole_leaf_index':sha(raw)},
                remaining=remaining, rows=records)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    for name in ('index_report', 'index', 'radial_certificate', 'emit'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--target', type=F, required=True)
    parser.add_argument('--radial-target', type=F, required=True)
    print(json.dumps(prepare(parser.parse_args()), separators=(',', ':')))
