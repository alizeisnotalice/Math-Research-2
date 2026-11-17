#!/usr/bin/env python3
"""Execute verify_switch_closure.py under Python 3.9 without changing its logic."""
import sys
from pathlib import Path

if len(sys.argv) < 2:
    raise SystemExit("usage: run_switch_closure_py39.py ORIGINAL [args ...]")
source = Path(sys.argv[1]).resolve()
text = source.read_text()
old = "(((x & 7) & row['permutation'][x >> 3]).bit_count() & 1)"
new = "(bin((x & 7) & row['permutation'][x >> 3]).count('1') & 1)"
if old not in text or "block.bit_count()" not in text:
    raise AssertionError("unexpected verifier source")
text = text.replace(old, new).replace("block.bit_count()", "bin(block).count('1')")
sys.argv = [str(source), *sys.argv[2:]]
scope = {"__name__": "__main__", "__file__": str(source)}
exec(compile(text, str(source), "exec"), scope)
