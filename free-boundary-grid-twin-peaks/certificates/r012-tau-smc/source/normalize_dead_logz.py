#!/usr/bin/env python3
"""Repair only the documented non-JSON spelling of a zero-weight dead path.

The old emitter wrote ``logZ:-inf``.  This program accepts that spelling only
when the entire record is an explicitly dead, zero-Z, zero-U record, changes
only that token to JSON ``null``, and rejects every other occurrence.
"""
import json
from pathlib import Path
import sys

src, dst = map(Path, sys.argv[1:3])
out = []
changed = 0
needle = '"logZ":-inf'
required = '"Z":0.0,"U":0.0,"min_distinct_ancestors":0,"dead":true'
for number, line in enumerate(src.read_text().splitlines(), 1):
    if needle in line:
        if line.count(needle) != 1 or required not in line:
            raise SystemExit(f"refuse line {number}: nonstandard logZ is not an explicit zero dead path")
        line = line.replace(needle, '"logZ":null')
        changed += 1
    json.loads(line)
    out.append(line)
dst.write_text("\n".join(out) + "\n")
print(f"converted={changed} lines={len(out)}")
