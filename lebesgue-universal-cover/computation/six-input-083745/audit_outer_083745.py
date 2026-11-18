"""Check the .83745 canonical seed and compare all fresh outer-screen packets."""
import argparse
from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
from io import BytesIO, StringIO
import json
from pathlib import Path
import subprocess
import sys
from tarfile import open as open_tar
from types import SimpleNamespace

import numpy as np
from flint import arb, ctx


TARGET = Fraction(16749, 20000)
COORD_DENOMINATOR = 2**44
AREA_DENOMINATOR = 2**57
GOAL = int(TARGET * AREA_DENOMINATOR) + 1
PARTITIONS = [431781, 115444, 17877, 6071, 3273]


def digest(raw):
    return sha256(raw).hexdigest()


def rational(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def exact_hull(points):
    rows = [(int(x), int(y), i) for i, (x, y) in enumerate(points)]
    rows.sort()
    rows = [row for i, row in enumerate(rows) if i == 0 or row[:2] != rows[i-1][:2]]

    def half(chain):
        result = []
        for p in chain:
            while len(result) >= 2:
                a, b = result[-2:]
                cross = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
                if cross > 0:
                    break
                result.pop()
            result.append(p)
        return result

    cycle = half(rows)[:-1] + half(reversed(rows))[:-1]
    return SimpleNamespace(vertices=np.array([row[2] for row in cycle], dtype=np.int64))


def replay_radial(comp):
    path = comp / "certify_radial_affine_fans.py"
    source = path.read_text()
    imported = "from scipy.spatial import ConvexHull"
    call = "ConvexHull(p.astype(float))"
    assert source.count(imported) == source.count(call) == 1
    source = source.replace(imported, "").replace(call, "ConvexHull(p)")
    with redirect_stdout(StringIO()) as output:
        exec(compile(source, str(path), "exec"), {"__file__": str(path), "ConvexHull": exact_hull})
    return json.loads(output.getvalue())


def source_tree_audit(comp, tree, seeds, index_dir, table32, table64):
    from certify_polygon_083731 import lower_radius
    from certify_refined_material_contraction import verify as verify_materials
    from select_tomar_083731_scores import radius_rows
    from tomar_d3_mask import CUT, leaf_mask

    seed_raw = (comp / "tomar-083745-future-seed.json").read_bytes()
    seed = json.loads(seed_raw)
    index_raw = (comp / "tomar-multi-target-after-083743.json").read_bytes()
    index = json.loads(index_raw)
    radial_raw = (comp / "five-body-radial-affine.json").read_bytes()
    radial = json.loads(radial_raw)
    masks = json.loads((comp / "five-orbit-slice.json").read_bytes())
    manifest_raw = (tree / "SHA256SUMS").read_bytes()
    pins = {name.removeprefix("./"): stamp
            for stamp, name in (line.split() for line in manifest_raw.decode().splitlines())}

    assert Fraction(seed["target"]) == Fraction(seed["radial_target"]) == TARGET
    assert seed["target_integer"] == GOAL
    assert seed["provenance"]["whole_leaf_index"] == digest(index_raw)
    assert seed["provenance"]["radial_certificate"] == digest(radial_raw)
    assert seed["tree_manifest_sha256"] == masks["manifest"] == digest(manifest_raw)
    assert index["provenance"]["manifest"] == digest(manifest_raw)
    assert index["area_denominator"] == AREA_DENOMINATOR
    assert index["coord_denominator"] == COORD_DENOMINATOR
    assert Fraction(index["max_index_target"]) == TARGET
    assert index["canonical_leaves"] == 6710314
    assert len(index["rows"]) == len(seed["rows"]) == 1024
    radial_replay = replay_radial(comp)
    assert radial_replay == radial
    signed_identity = json.loads(subprocess.check_output(
        [sys.executable, "-I", "-B", comp / "check_tomar_signed_translation_identity.py"]))
    materials = verify_materials(comp, tree, table32, table64)
    assert materials == json.loads((comp / "refined-material-contraction.json").read_bytes())

    root = json.loads((tree / "part0000.json").read_bytes())
    low, high = [list(map(Fraction, root[key])) for key in ("lo", "hi")]
    assert root["ns"] == [3, 5, 7, 9] and root["split0"] == [2] * 11
    assert all((v * COORD_DENOMINATOR).denominator == 1 for v in low + high)
    assert all((COORD_DENOMINATOR * (a + b) / 2).denominator == 1
               for a, b in zip(low, high, strict=True))
    assert all(low[axis] == 0 for axis in (2, 5, 8))
    pi = arb.pi()
    assert pi / 15 < rational(CUT) <= rational(Fraction(3684498829194, COORD_DENOMINATOR))
    assert rational(Fraction(3684498829194, COORD_DENOMINATOR)) < rational(high[2] / 2)
    assert all(rational(high[axis]) > 2*pi/n for axis, n in ((5, 7), (8, 9)))

    assert radius_rows(radial, TARGET) == seed["radii"]
    scale = Fraction(materials["curve_scale"])
    exterior = []
    for i, row in enumerate(seed["radii"]):
        k = 3 * i
        edge = min(-low[k], -low[k+1], high[k], high[k+1])
        assert Fraction(row["integer_radius"], COORD_DENOMINATOR) < Fraction(".069") < edge
        assert Fraction(row["integer_area_lower"], AREA_DENOMINATOR) >= Fraction(GOAL, AREA_DENOMINATOR)
        r = rational(scale / 2)
        other = rational(scale) * lower_radius(row["n"], pi)
        minimum = pi * (r**2 + other**2) / 2 + rational(edge) * (r + other)
        assert minimum > rational(Fraction(GOAL, AREA_DENOMINATOR))
        exterior.append(str(minimum))

    index_rows = {row[0]: row for row in index["rows"]}
    seed_rows = {row[0]: row for row in seed["rows"]}
    nodes = leaf_count = 0
    totals = Counter()
    active = set()
    for part, mask_info in enumerate(masks["rows"]):
        assert mask_info[0] == part
        stem = f"part{part:04}"
        head = (tree / (stem + ".json")).read_bytes()
        raw = (tree / (stem + ".npy")).read_bytes()
        assert digest(head) == pins[stem + ".json"] and digest(raw) == pins[stem + ".npy"]
        header = json.loads(head)
        assert header["part"] == part
        assert all(header[k] == root[k] for k in ("ns", "lo", "hi", "split0", "nparts"))
        code = np.load(BytesIO(raw), allow_pickle=False)
        assert code.dtype == np.uint8 and code.shape == (header["boxes"],)
        branch = code < 11
        assert np.all(branch | (code >= 253))
        need = 1 + np.cumsum(2 * branch.astype(np.int64) - 1)
        assert np.all(need[:-1] > 0) and need[-1] == 0
        nodes += len(code)
        nleaves = int(np.count_nonzero(~branch))
        leaf_count += nleaves
        if (part >> 8) & 1:
            assert part not in seed_rows and part not in index_rows
            expected_mask = np.where(branch, 0, 3).astype(np.uint8)
            assert mask_info[1:4] == [0, 0, nleaves]
            assert digest(expected_mask.tobytes()) == mask_info[-1]
            totals["outside"] += nleaves
            continue
        active.add(part)
        labels = leaf_mask(code, part, high[2])
        inside, crossing, outside = (int(np.count_nonzero(labels == k)) for k in (1, 2, 3))
        assert [inside, crossing, outside] == mask_info[1:4]
        assert digest(labels.tobytes()) == mask_info[-1]
        totals.update(inside=inside, crossing=crossing, outside=outside)
        with np.load(seeds / (stem + ".npz"), allow_pickle=False) as packet:
            positions, scores = packet["pos"], packet["score"]
        with np.load(index_dir / (stem + ".npz"), allow_pickle=False) as packet:
            reference_positions, reference_scores = packet["pos"], packet["full_floor"]
        assert positions.dtype == np.int32 and scores.dtype == np.int64
        assert np.array_equal(positions, np.unique(positions))
        assert np.all((labels[positions] == 1) | (labels[positions] == 2))
        assert np.all(scores < GOAL)
        assert np.array_equal(positions, reference_positions) and np.array_equal(scores, reference_scores)
        canonical, discarded, nonradial, remain = seed_rows[part][1:5]
        assert canonical == inside + crossing == index_rows[part][1]
        assert discarded == canonical - len(positions) and nonradial == 0 and remain == len(positions)
        assert [digest(positions.tobytes()), digest(scores.tobytes())] == seed_rows[part][-2:] == index_rows[part][-2:]
        totals["inherited_excluded"] += discarded
        totals["seed"] += remain
    assert active == set(seed_rows) == set(index_rows)
    assert totals["inside"] + totals["crossing"] == index["canonical_leaves"]
    assert nodes == 2 * leaf_count - 2048
    assert totals["seed"] == seed["remaining"] == PARTITIONS[0]
    assert totals["inherited_excluded"] == 6278533
    return {
        "partitions": 2048, "active_partitions": len(active), "original_nodes": nodes,
        "original_leaves": leaf_count, "mask_totals": dict(totals),
        "goal_integer": GOAL, "seed_sha256": digest(seed_raw), "index_sha256": digest(index_raw),
        "material_squared_maxima": materials["material_squared_maxima"],
        "radial_exterior_lower_intervals": exterior, "radial_replay": radial_replay["rows"],
        "signed_area_identity": signed_identity,
        "historic_index_geometry": "inherited .83743 whole-leaf scores and original index truncation",
    }


def compare_packets(package, seeds, replay):
    ledger = json.loads((replay / "stage-ledger.json").read_bytes())
    assert ledger["complete"] and ledger["parts"] == 1024 and ledger["integer_goal"] == GOAL
    kernel = replay / "outer-scored.cpp"
    kernel_sha = digest(kernel.read_bytes())
    assert kernel_sha == ledger["source_sha"] == "8ba4c735a230de2c395bf219ce674f2d4fdcdcdef4d07e1986f088a16d9469a6"
    summaries = []
    survivors = {}
    ought = {p for p in range(2048) if (p >> 8) & 1 == 0}
    for stage, stats in enumerate(ledger["audits"], 1):
        parts = dict(stats["work_files"])
        assert set(parts) == ought
        totals = Counter()
        with open_tar(package / f"native45-screen-{stage}.tar.gz", "r:gz") as bundle:
            inventory = dict((name, (size, stamp))
                             for name, size, stamp in json.load(bundle.extractfile("MANIFEST.json")))
            assert len(inventory) == 2048
            for member in bundle:
                if member.name == "MANIFEST.json":
                    continue
                size, stamp = inventory.pop(member.name)
                data = bundle.extractfile(member).read()
                assert len(data) == size and digest(data) == stamp
                path = replay / member.name.removeprefix("screen/")
                assert path.read_bytes() == data
                if path.suffix != ".json":
                    continue
                receipt = json.loads(data)
                config = receipt["config"]
                part = config["part"]
                assert config["layer"] == stage and config["kernel_sha256"] == kernel_sha
                assert config["goal"] == GOAL and Fraction(config["radial_goal"]) == TARGET
                assert receipt["data_sha256"] == parts[part]
                input_dir = seeds if stage == 1 else replay / f"stage{stage-1}"
                assert digest((input_dir / (path.stem + ".npz")).read_bytes()) == config["source_sha256"]
                with np.load(path.with_suffix(".npz"), allow_pickle=False) as packet:
                    pending = packet["remaining"]
                    accepted = packet["accepted"]
                    assert np.all(packet["full_floor"] >= GOAL)
                    if stage == 4 and len(pending):
                        survivors[str(part)] = pending.tolist()
                with np.load(input_dir / (path.stem + ".npz"), allow_pickle=False) as packet:
                    assert np.array_equal(np.sort(np.r_[accepted, pending]),
                                          packet["pos" if stage == 1 else "remaining"])
                for key in ("tested", "accepted", "remaining"):
                    totals[key] += receipt[key]
                totals["calls"] += receipt["counts"][0]
            assert not inventory
        assert all(stats[key] == value for key, value in totals.items())
        assert [totals["tested"], totals["remaining"]] == PARTITIONS[stage-1:stage+1]
        summaries.append(dict(totals))
    assert survivors == json.loads((replay / "retained.json").read_bytes())

    with open_tar(package / "native45-receipts.tar.gz", "r:gz") as bundle:
        scope = json.load(bundle.extractfile("inputs/ordinary/scope.json"))
    plain, deferred = (set(map(tuple, scope[key])) for key in ("leaves", "deferred"))
    final = {(int(part), pos) for part, positions in survivors.items() for pos in positions}
    assert not plain & deferred and plain | deferred == final
    return {"kernel_sha256": kernel_sha, "byte_identical_stage_packets": 8192,
            "fresh_stages": summaries, "fresh_retained_sha256": digest((replay / "retained.json").read_bytes()),
            "final_ordinary_receipt_leaves": len(plain), "final_deferred_receipt_leaves": len(deferred)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("seed", "screen"))
    for name in ("computation", "tree", "seeds", "whole_leaf_index", "table32",
                 "table64", "package", "screen"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    assert __debug__
    sys.path.insert(0, str(args.computation.resolve()))
    ctx.prec = 160
    if args.mode == "seed":
        result = source_tree_audit(args.computation, args.tree, args.seeds,
                                   args.whole_leaf_index, args.table32, args.table64)
    else:
        result = compare_packets(args.package, args.seeds, args.screen)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
