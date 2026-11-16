#!/usr/bin/env python3
"""Small exact regression for the four-block Franel / terminal canonicalization count."""
import importlib.util
from pathlib import Path

source = Path(__file__).with_name("verify_coprime_v2_generic.py")
spec = importlib.util.spec_from_file_location("generic_verify", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# F_4(2)=18, F_4(3)=164.  Nonterminal E=2,3,4,5 have radix products
# 3^4, 2^4, 3^4, 6^4, and terminal E=5 has canonical 2^(4-1)=8.
assert module.franel(2, 4) == 18
assert module.franel(3, 4) == 164
assert module.expected_layout_total(2, 3, 4) == 1482
assert 3**4 + 2**4 + 3**4 + 6**4 + 2**3 == 1482
assert module.expected_layout_total(2, 3, 5) == 8310
assert 3**5 + 2**5 + 3**5 + 6**5 + 2**4 == 8310
assert module.decode(17, [3, 3, 2]) == [2, 2, 1]
print("PASS: generic n=4/n=5 Franel formula and terminal 2^(n-1)-2 correction")
