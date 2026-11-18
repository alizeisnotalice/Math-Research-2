"""Bind stage-four retained keys to complete terminal receipts and replays."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
import tarfile

import numpy as np


def digest(raw):
    return sha256(raw).hexdigest()


def bundle(archive):
    records = {}
    with tarfile.open(archive, "r:gz") as packed:
        manifest = json.load(packed.extractfile("MANIFEST.json"))
        assert len(manifest) + 1 == len(set(packed.getnames())) == len(packed.getnames())
        for filename, length, checksum in manifest:
            raw = packed.extractfile(filename).read()
            assert len(raw) == length and digest(raw) == checksum and filename not in records
            records[filename] = raw
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("package", "comp", "tree", "dense_replay", "ordinary_matched", "output"):
        parser.add_argument(field, type=Path)
    parser.add_argument("--pure-replay", type=Path, required=True)
    args = parser.parse_args()
    sys.path[:0] = [str(args.package.parent), str(args.comp)]
    from certify_tomar_six_high import FINE, SourceTree
    from verify_five_selected_083740_traces import audit

    original = SourceTree(args.tree)
    goal = (Fraction(".83745") * 2**57).__ceil__()
    assert goal == 120689264254125501
    reference = bundle(args.package / "native45-receipts.tar.gz")
    retained = bundle(args.package / "native45-screen-4.tar.gz")
    stage_index_raw = reference["inputs/ordinary/stage4-parts.json"]
    stage_rows = json.loads(stage_index_raw)
    stage_index = dict(stage_rows)
    assert len(stage_rows) == len(stage_index) == 1024
    scope = json.loads(reference["inputs/ordinary/scope.json"])
    assert digest(stage_index_raw) == scope["parent_parts_sha256"]
    ordinary_scope = set(map(tuple, scope["leaves"]))
    deferred_scope = set(map(tuple, scope["deferred"]))
    assert len(ordinary_scope) == len(scope["leaves"]) == 3071
    assert len(deferred_scope) == len(scope["deferred"]) == 202
    assert not ordinary_scope & deferred_scope
    stage_keys = set()
    for part, index_hash in stage_index.items():
        prefix = f"screen/stage4/part{part:04d}"
        raw = retained[prefix + ".json"]
        header = json.loads(raw)
        data = retained[prefix + ".npz"]
        assert digest(raw) == index_hash
        assert digest(data) == header["data_sha256"]
        config = header["config"]
        assert config["part"] == part and config["layer"] == 4 and config["goal"] == goal
        with np.load(BytesIO(data), allow_pickle=False) as npz:
            keys = [(part, int(pos)) for pos in npz["remaining"]]
        assert len(keys) == len(set(keys)) == header["remaining"]
        stage_keys.update(keys)
    assert len(stage_keys) == 3273 and stage_keys == ordinary_scope | deferred_scope

    replayed = {}
    for part, pos, filename, checksum in json.loads(args.ordinary_matched.read_bytes()):
        assert (part, pos) not in replayed
        raw = Path(filename).read_bytes()
        assert digest(raw) == checksum
        replayed[part, pos] = raw

    terminal = {}
    ordinary_failed = set()
    ordinary_n = 0
    ordinary_nodes = 0
    ordinary_matched = set()
    for name, raw in sorted(reference.items()):
        if not name.startswith("ordinary/") or "/rows/" not in name:
            continue
        row = json.loads(raw)
        key = row["part"], row["pos"]
        ordinary_n += 1
        stat = row["stat"]
        ordinary_nodes += stat[1]
        if key in replayed:
            assert replayed[key] == raw
            ordinary_matched.add(key)
        assert key in ordinary_scope and key not in terminal and key not in ordinary_failed
        box = original.leaf(*key)
        assert all(a < b and (a * FINE).denominator == (b * FINE).denominator == 1
                   for a, b in box)
        if stat[0]:
            assert stat[5] == 0 and stat[1] == 2*(stat[2] + stat[3]) - 1
            terminal[key] = ["ordinary", stat[1], digest(raw), key in ordinary_matched]
        else:
            assert stat[5] == 2 and stat[4] == row["depth"] == 135
            ordinary_failed.add(key)
    assert ordinary_n == 3071 and len(terminal) == 3053 and len(ordinary_failed) == 18
    assert set(replayed) == ordinary_matched

    pure_raw = args.pure_replay.read_bytes()
    pure_replay = json.loads(pure_raw)
    assert pure_replay["complete"] and Fraction(pure_replay["threshold"]) == Fraction(".837455")
    assert pure_replay["original_leaves"] == pure_replay["original_record_matches"] == pure_replay["original_trace_matches"] == 196
    assert pure_replay["replayed_geometric_visits"] == 348400
    pure_index = {(part, pos): (visits, ends, same_result, same_trace)
                  for part, pos, visits, ends, same_result, same_trace in pure_replay["rows"]}
    assert len(pure_replay["rows"]) == len(pure_index) == 196
    checked_pure = set()
    restored = {}
    for volume in ("native45-pure-a.tar.gz", "native45-pure-b.tar.gz"):
        data = bundle(args.package / volume)
        assert not restored.keys() & data.keys()
        restored.update(data)
    for name, raw in sorted(restored.items()):
        if not name.startswith("pure/actual/0.837455/") or not name.endswith(".json"):
            continue
        row = json.loads(raw)
        key = row["part"], row["pos"]
        assert key in stage_keys and key not in terminal
        assert key not in ordinary_failed
        assert row["original_header_sha256"] == original.pins[f"part{key[0]:04}.json"]
        assert row["original_tree_sha256"] == original.pins[f"part{key[0]:04}.npy"]
        assert Fraction(row["threshold"]) == Fraction(".837455")
        box = original.leaf(*key)
        root = [int(value * FINE) for vector in zip(*box) for value in vector]
        assert root == row["root"]
        trace = restored[name[:-5] + ".bin"]
        assert row["trace_sha256"] == digest(trace)
        counters = audit(root, trace)
        assert counters == {key: row["result"][key] for key in counters}
        assert not row["result"]["remaining"]
        visits, ends, identical_result, identical_trace = pure_index[key]
        assert visits == counters["vis"] and ends == sum(counters["closed"].values())
        assert identical_result and identical_trace
        checked_pure.add(key)
        terminal[key] = ["pure_five", counters["vis"], digest(trace), True]
    assert checked_pure == pure_index.keys()

    dense_keys = set()
    for volumes, prefix, kind in ((reference, "dense/native45-dense-", "dense_rescue"),
                                  (restored, "dense/proof-", "dense_original")):
        for name, raw in sorted(volumes.items()):
            if not name.startswith(prefix) or not name.endswith("/result.json"):
                continue
            row = json.loads(raw)
            key = row["part"], row["pos"]
            assert key in stage_keys and key not in terminal
            assert (key in ordinary_failed) == (kind == "dense_rescue")
            assert row["all_closed"] and row["complete"] and row["threshold_integer"] == goal
            assert Fraction(row["target"]) == Fraction(".83745") and len(row["rows"]) == 1
            leaf = row["rows"][0]
            lo, hi = [[int(v * FINE) for v in vector] for vector in zip(*original.leaf(*key))]
            assert (lo, hi) == (row["original_lo"], row["original_hi"]) == (leaf["lo"], leaf["hi"])
            assert leaf["bits"] == "" and leaf["terminal_guard"] == "mesh64-dense"
            stat = leaf["stat"]
            assert stat[0] == 1 and stat[5] == 0 and stat[1] == 2*(stat[2] + stat[3]) - 1
            second_pass = (args.dense_replay / f"{key[0]}-{key[1]}" / "result.json").read_bytes()
            assert json.loads(second_pass) == row
            terminal[key] = [kind, stat[1], digest(raw), True]
            dense_keys.add(key)
    assert len(dense_keys) == 24 and dense_keys >= ordinary_failed
    assert terminal.keys() == stage_keys

    counts, visits, replay_counts = Counter(), Counter(), Counter()
    for kind, nodes, _, checked in terminal.values():
        counts[kind] += 1
        visits[kind] += nodes
        replay_counts[kind] += checked
    assert counts == {"ordinary": 3053, "dense_rescue": 18, "dense_original": 6, "pure_five": 196}
    report = {
        "target": ".83745",
        "stage_four_count": len(stage_keys),
        "terminal_candidate_counts": counts,
        "terminal_node_counts": visits,
        "terminal_nodes_total": sum(visits.values()),
        "pure_five_replay_manifest_sha256": digest(pure_raw),
        "terminal_replay_counts": replay_counts,
        "all_terminal_replayed": counts == replay_counts,
        "all_original_ordinary_replayed": ordinary_n == len(ordinary_matched),
        "ordinary_including_failed": {"first_pass": ordinary_n, "nodes": ordinary_nodes,
                                     "second_pass": len(ordinary_matched),
                                     "matched_failed": sorted(ordinary_failed & ordinary_matched)},
        "leaf_fields": ["part", "pos", "terminal_method", "nodes", "original_certificate_sha256",
                        "second_pass_verified"],
        "leaves": [[*leaf, *terminal[leaf]] for leaf in sorted(stage_keys)],
    }
    args.output.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "leaves"}, indent=2))


if __name__ == "__main__":
    assert __debug__
    main()
