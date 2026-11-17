#!/usr/bin/env python3
"""Read-only audit of the archived (L,d,e) tables at L=2,3.

No DP is rerun.  H[L,d,e] is the punctured/stable-fibre count; the script
tests both the unweighted global identity H_2=2H_3 and its proposed
distance/squared-distance refinements.
"""

from collections import defaultdict
import hashlib
import json
from pathlib import Path


import argparse
_parser=argparse.ArgumentParser()
_parser.add_argument('--geometry-dir',type=Path,required=True)
_args=_parser.parse_args()
SOURCES={3:_args.geometry_dir/'gap-n3-all.json',5:_args.geometry_dir/'n5-gap-geometry-summary.json'}


def aggregate(source: Path) -> dict[tuple[int, int, int], int]:
    data = json.loads(source.read_text())
    if 'holes' in data:
        rows = [row for hole in data['holes'] for row in hole['joint_H_L_d_e2']]
    else:
        rows = data['joint_H_L_d_e2']
    ans: dict[tuple[int, int, int], int] = defaultdict(int)
    for L, d, e, count in rows:
        ans[int(L), int(d), int(e)] += int(count)
    return ans


def audit_one(n: int, source: Path) -> dict:
    h = aggregate(source)
    keys = sorted({(d, e) for L, d, e in h if L in (2, 3)})
    delta = {(d, e): h[2, d, e] - 2 * h[3, d, e] for d, e in keys}
    assert sum(delta.values()) == 0, (n, 'global h2 identity')
    assert all(d <= 4 for (d, e), z in delta.items() if z), (n, 'far support')
    nonzero = [
        {'d': d, 'e2': e, 'h2': str(h[2, d, e]), 'h3': str(h[3, d, e]),
         'h2_minus_2h3': str(delta[d, e])}
        for d, e in keys if delta[d, e]
    ]
    by_d: dict[int, int] = defaultdict(int)
    for (d, _e), z in delta.items():
        by_d[d] += z
    return {
        'n': n,
        'source': source.name,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'h2_total': str(sum(h[2, d, e] for d, e in keys)),
        'h3_total': str(sum(h[3, d, e] for d, e in keys)),
        'global_h2_minus_2h3': str(sum(delta.values())),
        'squared_distance_weighted_delta': str(sum(e * z for (_d, e), z in delta.items())),
        'manhattan_distance_weighted_delta': str(sum(d * z for (d, _e), z in delta.items())),
        'by_distance_delta': {str(d): str(z) for d, z in sorted(by_d.items()) if z},
        'nonzero_joint_cells': nonzero,
    }


def main() -> None:
    result = {'status': 'passed', 'scope': 'read-only exact finite audit; no asymptotic inference',
              'results': [audit_one(n, source) for n, source in SOURCES.items()]}
    out = Path(__file__).with_name('joint_h2_h3_audit.json')
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
