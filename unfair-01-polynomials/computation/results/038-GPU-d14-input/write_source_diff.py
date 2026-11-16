#!/usr/bin/env python3
"""Write an exact, machine-produced d13-to-d14 worker unified diff."""
import difflib
from pathlib import Path
HERE=Path(__file__).resolve().parent
OLD=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/ds_optimize/round033-gpu-d13/worker/gpu_factor_search_prefix.py')
NEW=HERE/'worker/gpu_factor_search_prefix.py'
out=''.join(difflib.unified_diff(OLD.read_text().splitlines(True),NEW.read_text().splitlines(True),fromfile=str(OLD),tofile=str(NEW)))
assert out and 'd <= 13' in out and 'd <= 14' in out
(HERE/'source.diff').write_text(out)
print(len(out.encode()))
