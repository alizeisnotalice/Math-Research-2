"""Restore the exact historical reference at its recorded temp path if absent.
Never replace an existing file with different bytes.
"""
from pathlib import Path
import hashlib,os
assert os.path.ismount('/Volumes/data')
r=Path(__file__).resolve().parent
source=r.parent/'038-GPU-d13/worker/gpu_factor_search_prefix.py'
target=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/ds_optimize/round033-gpu-d13/worker/gpu_factor_search_prefix.py')
raw=source.read_bytes();expected='4cbe9ed5f04acebfaf6f420ce5a7c076b086352c317376cd14ef0520a8acff70'
assert hashlib.sha256(raw).hexdigest()==expected
if target.exists():assert target.read_bytes()==raw
else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
print('frozen source available with exact recorded hash',expected)
