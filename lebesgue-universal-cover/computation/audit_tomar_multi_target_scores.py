"""Build one authenticated whole-leaf integer index for a range of area targets."""
from argparse import ArgumentParser
import ctypes
from fractions import Fraction as F
from io import BytesIO
import json
from math import isqrt
import os
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

from flint import arb, ctx
import numpy as np

from audit_tomar_outer_083730 import AREA, FINE, ceil, score_floor, sha
from select_tomar_083731_scores import q, radius_rows

# Exact full-domain integer from proofs/tomar-five-083729.md.
GLOBAL_FLOOR = 120666205824033364
TARGETS = ('.83731', '.83732', '.83733', '.83734', '.83735', '.83736', '.83740', '.83745')
GAP_EDGES = tuple(map(F, ('.0000001', '.000001', '.000002', '.000005', '.00001',
                          '.00002', '.00003', '.00004', '.00006', '.00010', '.00017')))


def checked_json(path, expected=None):
    payload = path.read_bytes()
    digest = sha(payload)
    if expected:
        assert digest == expected
    return json.loads(payload), digest


def library(work, cache):
    source = Path(__file__).with_name('tomar_multi_target_radial.cpp')
    binary = work/'tomar-multi-target-radial.dylib'
    subprocess.run([os.environ.get('CXX', 'clang++'), '-O3', '-std=c++17', '-shared', '-fPIC',
                    '-ftrapv', f'-fmodules-cache-path={cache}', str(source), '-o', str(binary)], check=True)
    handle = ctypes.CDLL(str(binary))
    f = handle.radial_scores
    i64 = np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags='C_CONTIGUOUS')
    u8 = np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1, flags='C_CONTIGUOUS')
    f.argtypes = (u8, ctypes.c_int, i64, i64, u8, i64, i64, i64)
    return f, sha(source.read_bytes())


def radial_floor(distance, params):
    rn, rd, low, up = map(int, params[:4])
    excess = distance*rd-rn*FINE
    return low + excess*(low-up)//(rn*FINE) if excess >= 0 else 0


def support(target, gap):
    # The target minus the leaf floor is strictly less than gap/AREA.
    required = ceil(4*target*F(gap, AREA)*10**16)
    s = isqrt(required)
    return s + int(s*s < required)


def main(a):
    ctx.prec = 160
    sys.path.insert(0, str(a.computation))
    from tomar_d3_mask import CUT, leaf_mask

    a.work.mkdir(parents=True, exist_ok=True)
    scanner, source_digest = library(a.work, a.cache)
    files = {}
    for key, filename in {
        'bucket_index': 'tomar-grouped-08371-index.json',
        'bucket_ledger': 'tomar-grouped-08371-ledger.json',
        'canonical_slice': 'five-orbit-slice.json',
        'prior_radial': 'five-body-radial-affine.json',
        'stage_index': 'tomar-083731-score-archive-index.json',
    }.items():
        files[key] = checked_json(a.computation/filename)
    (bidx, bsha), (bled, blsha), (slices, ssha), (oldrad, rsha), (lidx, lsha) = files.values()
    assert bidx['audit_sha256'] == blsha and bled['slice_index_sha256'] == ssha
    assert bidx['complete'] and bled['complete'] and not bled['historical_area_inputs']
    assert bidx['target'] == bled['target'] == '0.8371'
    assert oldrad['source_sha256'] == sha((a.computation/'certify_radial_affine_fans.py').read_bytes())
    assert oldrad['table_sha'] == sha((a.computation/'five-body-radial-material.npz').read_bytes())
    stage, stagehash = checked_json(a.stages/'stage-ledger.json', lidx['stage_ledger_sha256'])
    old_tail, tailsha = checked_json(a.stages/'retained.json', lidx['retained_sha256'])
    cuts = [{part: digest for part, digest in row['work_files']} for row in stage['audits']]
    bucket_receipts = {p: pair[:2] for p, *pair in bled['receipt_hashes']}
    assert stage['complete'] and stage['parts'] == len(bucket_receipts) == 1024
    assert len(cuts) == 4 and all(set(c) == set(bucket_receipts) for c in cuts)
    assert sum(map(len, old_tail.values())) == 614

    manifest = (a.tree/'SHA256SUMS').read_bytes()
    manifest_sha = sha(manifest)
    assert manifest_sha == bled['input_tree_sha256'] == slices['manifest']
    tree_files = {name.removeprefix('./'): digest
                  for digest, name in (row.split() for row in manifest.decode().splitlines())}
    root, rootsha = checked_json(a.tree/'part0000.json', tree_files['part0000.json'])
    lo, hi = (np.array([int(F(v)*FINE) for v in root[k]], np.int64) for k in ('lo', 'hi'))
    assert all(F(int(i), FINE) == F(v) for k, arr in (('lo', lo), ('hi', hi))
               for i, v in zip(arr, root[k], strict=True))
    assert np.all((lo+hi) % 2 == 0) and np.all(lo < hi)
    mid = (lo+hi)//2
    assert root['split0'] == [2]*11 and root['ns'] == [3, 5, 7, 9]
    assert all(lo[j] == 0 for j in (2, 5, 8))
    cxx_cut = F(3684498829194, FINE)
    assert arb.pi()/15 < q(CUT) <= q(cxx_cut) < q(F(int(hi[2]), FINE))
    assert all(q(F(int(hi[j]), FINE)) > 2*arb.pi()/n for j, n in ((5, 7), (8, 9)))

    queries = sorted(set(a.targets or map(F, TARGETS)))
    assert queries[0] >= F(stage['target'])
    goals = np.array([int(t*AREA)+1 for t in queries], np.int64)
    max_goal = goals[-1]
    radii = radius_rows(oldrad, queries[-1])
    old_clip = radius_rows(oldrad, F(stage['radial_target']))
    param_rows = []
    axis_edges = []
    for b, row in enumerate(radii):
        r = F(row['previous_radius'])
        param = (r.numerator, r.denominator, row['old_floor'], ceil(AREA*F(row['center_upper'])),
                 old_clip[b]['integer_radius'])
        assert radial_floor(param[-1], param) >= old_clip[b]['integer_area_lower']-1
        param_rows.append(param)
        edge = min(-lo[3*b], hi[3*b], -lo[3*b+1], hi[3*b+1])
        assert row['integer_radius'] < edge < FINE
        assert radial_floor(row['integer_radius'], param) >= max_goal
        score = radial_floor(int(edge), param)
        assert score >= max_goal
        axis_edges.append([3+2*b, int(edge), score])
    param = np.array(param_rows, np.int64).ravel()

    table = np.array([max(GLOBAL_FLOOR, score_floor(b)) if b else GLOBAL_FLOOR
                      for b in range(65536)], np.int64)
    stage_goal = stage['integer_goal']
    previous = np.zeros((len(queries), 6), np.int64)
    all_scores, partition_rows, partition_gaps = [], [], []
    seen = set()
    if a.emit:
        a.emit.mkdir(parents=True, exist_ok=True)
    for record in bidx['archives']:
        archive = a.computation/record['name']
        assert sha(archive.read_bytes()) == record['sha256']
        with ZipFile(archive) as z:
            for name in z.namelist():
                if not (name.startswith('native/') and name.endswith('.npz')):
                    continue
                stem = Path(name).stem
                part = int(stem[4:])
                assert part not in seen
                seen.add(part)
                packet, receipt = z.read(name), z.read(name.replace('.npz', '.json'))
                assert [sha(receipt), sha(packet)] == bucket_receipts[part]
                meta = json.loads(receipt)
                code_bytes = (a.tree/(stem+'.npy')).read_bytes()
                header, head_sha = checked_json(a.tree/(stem+'.json'), tree_files[stem+'.json'])
                assert head_sha == meta['header_sha256']
                assert sha(code_bytes) == tree_files[stem+'.npy'] == meta['codes_sha256']
                assert all(root[k] == header[k] for k in ('lo', 'hi', 'ns', 'split0'))
                codes = np.load(BytesIO(code_bytes), allow_pickle=False)
                assert codes.dtype == np.uint8 and codes.ndim == 1
                mask = (np.load(a.masks/(stem+'.npy'), allow_pickle=False) if a.masks else
                        leaf_mask(codes, part, F(root['hi'][2])))
                assert mask.dtype == np.uint8 and mask.shape == codes.shape
                assert sha(mask.tobytes()) == slices['rows'][part][-1] == meta['slice_mask_sha256']
                active = ((mask == 1)|(mask == 2)).astype(np.uint8)
                pos = np.flatnonzero(active).astype(np.int32)
                with np.load(BytesIO(packet), allow_pickle=False) as data:
                    buckets = data['scores']
                assert buckets.dtype == np.uint16 and len(buckets) == len(pos) == meta['all_original']

                pre_radial = table[buckets]
                combined = pre_radial.copy()
                inside = np.zeros(len(pos), np.int64)
                for j, part_hash in enumerate(cuts, start=1):
                    packet = (a.stages/f'stage{j}'/(stem+'.npz')).read_bytes()
                    assert sha(packet) == part_hash[part]
                    with np.load(BytesIO(packet), allow_pickle=False) as data:
                        accepted, full, clipped = data['accepted'], data['full_floor'], data['clipped_floor']
                    assert accepted.dtype == np.int32 and full.dtype == np.int64 and accepted.shape == full.shape
                    idx = np.searchsorted(pos, accepted)
                    assert np.all(idx < len(pos)) and np.array_equal(pos[idx], accepted)
                    assert np.all(full >= stage_goal) and len(idx) == len(np.unique(idx))
                    assert clipped.dtype == np.int64 and clipped.shape == full.shape
                    assert np.all((clipped == 0)|(clipped >= full)) and not inside[idx].any()
                    combined[idx] = np.maximum(combined[idx], full)
                    inside[idx] = clipped
                bits = np.array(np.unravel_index(part, (2,)*11), bool)
                low, high = np.where(bits, mid, lo), np.where(bits, hi, mid)
                radial = np.zeros(len(pos), np.int64)
                exterior = np.zeros_like(radial)
                assert scanner(codes, len(codes), low, high, active, param, radial, exterior) == 0
                assert np.all(radial <= max(row[2] for row in axis_edges)+AREA)
                original_bounds = np.maximum(combined, radial)
                lifted = np.minimum(inside, exterior)
                bounds = np.maximum(original_bounds, lifted)
                rare = bounds < max_goal
                values = bounds[rare]
                residue = pos[rare]
                tail = np.isin(pos, np.array(old_tail.get(str(part), []), np.int32), assume_unique=True)
                assert np.all((bounds >= stage_goal)|tail)

                if a.emit:
                    np.savez_compressed(a.emit/(stem+'.npz'), pos=residue, full_floor=values)
                all_scores.append(values)
                counts = [int((values < goal).sum()) for goal in goals]
                partition_rows.append([part, len(pos), *counts, sha(residue.tobytes()), sha(values.tobytes())])
                for j, goal in enumerate(goals):
                    bucket_covered = pre_radial >= goal
                    affine_covered = combined >= goal
                    pre_lift_covered = original_bounds >= goal
                    final_covered = bounds >= goal
                    previous[j] += (bucket_covered.sum(), (affine_covered & ~bucket_covered).sum(),
                                    (pre_lift_covered & ~affine_covered).sum(), (final_covered & ~pre_lift_covered).sum(),
                                    (~final_covered).sum(),
                                    (tail & ~final_covered).sum())
                focus = int(F('.83735')*AREA)+1
                selected = values[values < focus]
                if selected.size:
                    partition_gaps.append([part, len(selected), int(focus-selected.min())])

    assert seen == set(bucket_receipts)
    total = sum(row[1] for row in partition_rows)
    assert total == 6710314 and np.all(previous[:, :5].sum(axis=1) == total)
    scores = np.sort(np.concatenate(all_scores))
    curve = []
    edges = [int(v*AREA) for v in GAP_EDGES]
    for t, gate, counts in zip(queries, goals, previous, strict=True):
        goal = int(gate)
        unresolved = int(np.searchsorted(scores, goal))
        assert unresolved == counts[4]
        gaps = goal-scores[:unresolved]
        # Counts use successive inclusive upper gap edges.
        cumulative = [unresolved-int(np.searchsorted(scores, goal-e, side='left')) for e in edges]
        histogram = np.diff([0, *cumulative, unresolved]).tolist()
        worst = int(gaps[0]) if unresolved else 0
        curve.append(dict(target=str(t), goal=goal, counts=counts.tolist(), gap_counts=histogram,
                          worst_gap_int=worst, radius_surplus_grid_1e8=support(t, worst)))
    if a.emit:
        np.save(a.emit/'sorted-floor.npy', scores, allow_pickle=False)
    focus_goal = int(F('.83735')*AREA)+1
    unique, frequency = np.unique(scores[scores < focus_goal], return_counts=True)
    groups = sorted(zip(map(int, frequency), map(int, unique), strict=True), reverse=True)[:16]
    return dict(area_denominator=AREA, coord_denominator=FINE, max_index_target=str(queries[-1]),
                canonical_leaves=total, prior_global_integer=GLOBAL_FLOOR,
                count_labels=['old_bucket', 'new_stage', 'new_radial', 'radial_recombined', 'unresolved', 'of_those_old_614'],
                radial_params=param_rows, radial_coordinates_for_max=radii, radial_root_edges=axis_edges,
                old_clip_rows=old_clip, angular_cxx_cut=str(cxx_cut),
                radial_coordinates_for_083735=radius_rows(oldrad, F('.83735')),
                gap_edges=[str(e) for e in GAP_EDGES], curve=curve,
                focus_largest_partitions=sorted(partition_gaps, key=lambda x: (-x[1], x[0]))[:16],
                focus_largest_integer_plateaus=groups,
                sorted_floor_sha256=sha(scores.tobytes()), rows=sorted(partition_rows),
                provenance={**{k:v[1] for k, v in files.items()}, 'manifest':manifest_sha, 'root':rootsha,
                            'stage_ledger':stagehash, 'stage_retained':tailsha, 'cpp':source_digest,
                            'python':sha(Path(__file__).read_bytes())})


def query(a):
    report = json.loads(a.report.read_bytes())
    max_target = F(report['max_index_target'])
    assert report['area_denominator'] == AREA
    if a.part is None:
        scores = np.load(a.index/'sorted-floor.npy', allow_pickle=False)
        assert sha(scores.tobytes()) == report['sorted_floor_sha256']
        assert scores.dtype == np.int64 and np.all(scores[:-1] <= scores[1:])
    else:
        row = next(row for row in report['rows'] if row[0] == a.part)
        with np.load(a.index/f'part{a.part:04}.npz', allow_pickle=False) as data:
            positions, scores = data['pos'], data['full_floor']
        assert positions.dtype == np.int32 and scores.dtype == np.int64 and scores.shape == positions.shape
        assert [sha(positions.tobytes()), sha(scores.tobytes())] == row[-2:]
    output = []
    for t in a.targets:
        if not 0 < t <= max_target:
            raise ValueError(f'target {t} lies outside the certified query range (0,{max_target}]')
        g = int(t*AREA)+1
        if a.part is None:
            size = int(np.searchsorted(scores, g))
            row = {'remaining': size}
            gaps = g-scores[:size]
        else:
            keep = scores < g
            row = {'part': a.part, 'remaining': int(keep.sum()),
                   'positions': positions[keep].tolist(), 'full_floor': scores[keep].tolist()}
            gaps = g-scores[keep]
        worst = int(gaps.max(initial=0))
        row.update(target=str(t), goal=g, worst_gap_int=worst, radius_surplus_grid_1e8=support(t, worst))
        output.append(row)
    return output


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'query':
        parser = ArgumentParser(description=query.__name__)
        parser.add_argument('report', type=Path)
        parser.add_argument('index', type=Path)
        parser.add_argument('targets', nargs='+', type=F)
        parser.add_argument('--part', type=int)
        print(json.dumps(query(parser.parse_args(sys.argv[2:])), separators=(',', ':')))
        sys.exit(0)

    parser = ArgumentParser(description=__doc__)
    for key in ('computation', 'tree', 'stages', 'work', 'cache'):
        parser.add_argument(key, type=Path)
    parser.add_argument('--masks', type=Path)
    parser.add_argument('--emit', type=Path)
    parser.add_argument('--targets', nargs='+', type=F)
    print(json.dumps(main(parser.parse_args()), separators=(',', ':')))
