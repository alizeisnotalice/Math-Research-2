#!/usr/bin/env python3
"""Independently authenticate completed degree-64 checkpoints and their affine residues."""

import argparse
import json
from pathlib import Path
import sys

from replay import P, ROOT, encoded, file_digest, verify_result


def main():
    if sys.flags.optimize:
        raise SystemExit("Python assertions must be enabled.")
    p=argparse.ArgumentParser()
    p.add_argument("run",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--compare-v1",type=Path)
    p.add_argument("--compare-deterministic",type=Path)
    p.add_argument("--export-deterministic",type=Path)
    args=p.parse_args()
    assert not args.out.exists()
    run=args.run
    manifest=json.loads((run/"manifest.json").read_bytes())
    from hashlib import sha256
    assert manifest["q"]==list(range(1,P))
    assert manifest["source_sha256"]==file_digest(ROOT/"ordinary_remainder.cpp")
    assert manifest["reference_source_sha256"]==file_digest(ROOT/"reference_source.cpp")
    assert manifest["runner_sha256"]==file_digest(ROOT/"replay.py")
    summary=json.loads((run/"complete.json").read_bytes())
    assert summary["parameters"]==162 and summary["free_bits"]==41
    assert summary["total_stored"]==162*(1<<20)
    assert summary["total_scanned"]==162*(1<<21)
    assert summary["total_cross_coordinate_columns"]==162*65
    assert summary["manifest_sha256"]==file_digest(run/"manifest.json")
    assert summary["sorted_q_index_sha256"]==file_digest(run/"all-q.jsonl")
    assert len(list((run/"by-q").glob("*")))==162
    rows=list(map(json.loads,(run/"all-q.jsonl").read_bytes().splitlines()))
    assert [row["q"] for row in rows]==manifest["q"]
    deterministic=[]
    expected_projection_pairs=0
    for row in rows:
        q=row["q"]
        checkpoint=run/"by-q"/f"q{q:03d}.json"
        assert row["checkpoint_sha256"]==file_digest(checkpoint)
        c=json.loads(checkpoint.read_bytes())
        assert c["q"]==q and c["exit_code"]==0
        report=c["report"]
        check=verify_result(report,64,q,"tight")
        assert c["verification"]==check and all(row[k]==v for k,v in check.items())
        for k in ("stored_patterns","scanned_patterns","projected_pairs","matches"):
            assert row[k]==report[k]
        if args.compare_v1:
            old=json.loads((args.compare_v1/"by-q"/checkpoint.name).read_bytes())["report"]
            assert all(report[k]==v for k,v in old.items())
        assert not row["matches"] and not check["max_square_degree7_char0_squarefree"]
        expected_projection_pairs+=report["projected_pairs"]
        deterministic.append({k:v for k,v in row.items() if k!="checkpoint_sha256"})
    assert summary["full_matches"]==0 and summary["full_scope_matches"]==0
    assert summary["exact_projection_pairs"]==expected_projection_pairs
    raw=b"".join(map(encoded,deterministic))
    if args.compare_deterministic:
        assert raw==args.compare_deterministic.read_bytes()
    if args.export_deterministic:
        with args.export_deterministic.open("xb") as file:
            file.write(raw)
    verification=dict(status="verified",parameters=len(rows),
                      projected_masks=expected_projection_pairs,
                      reciprocal_column_comparisons=summary["total_cross_coordinate_columns"],
                      deterministic_index_sha256=sha256(raw).hexdigest(),
                      finished_checks_digest=file_digest(run/"all-q.jsonl"),
                      compared_v1=bool(args.compare_v1))
    args.out.write_bytes(encoded(verification))
    print(json.dumps(verification,sort_keys=True))


if __name__=="__main__":
    main()
