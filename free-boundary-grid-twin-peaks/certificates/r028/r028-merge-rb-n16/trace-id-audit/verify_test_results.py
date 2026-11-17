#!/usr/bin/env python3
"""Validate stable semantics of a trace-ID test replay, ignoring child PIDs."""
import argparse
import json
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("result", type=Path)
args = ap.parse_args()
got = json.loads(args.result.read_text())
assert got["status"] == "passed"
assert got["scope"] == "trace ID validation only; production kernel unchanged"
tests = got["tests"]
assert len(tests) == 4
for record, ids, accepted in zip(
    tests, ((0, 0), (1, 1), (0, 1), (1, 0)), (True, False, False, False)
):
    assert tuple(record["ids"]) == ids
    assert record["expected_accepted"] is accepted
    assert (record["returncode"] == 0) is accepted
    assert record["exit_confirmed"] is True
print("trace-ID replay semantics verified")
