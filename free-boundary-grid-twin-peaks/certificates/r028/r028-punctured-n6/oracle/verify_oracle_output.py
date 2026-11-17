#!/usr/bin/env python3
"""Compare the deterministic mathematical fields of two oracle outputs."""
from __future__ import annotations

import json
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("usage: verify_oracle_output.py expected.json fresh.json")
expected, fresh = (json.loads(Path(path).read_text()) for path in sys.argv[1:])
for key in ("status", "method", "scope", "runs", "reference_files", "source_sha256"):
    assert fresh[key] == expected[key], key
assert fresh["status"] == "passed"
print(json.dumps({"passed": True, "checked_fields": ["status", "method", "scope", "runs", "reference_files", "source_sha256"]}, indent=2))
