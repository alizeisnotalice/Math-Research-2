"""Recover all outstanding original positions from a complete list of stronger-score receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
for key in ("parent","primary_ledger","native","output"):p.add_argument(key,type=Path)
a=p.parse_args()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
parent=json.loads(a.parent.read_bytes())
report=json.loads(a.primary_ledger.read_bytes())
assert report["target"]=="0.8372" and report["authenticated_parent_ledger_sha256"]==sha(a.parent)
assert report["parts"]==parent["parts"] and report["totals"]["canonical"]==parent["totals"]["canonical_leaves"]
assert [x[0] for x in report["higher_receipt_hashes"]]==[x[0] for x in parent["receipt_hashes"]]
left={}
for part,h1,h2 in report["higher_receipt_hashes"]:
    stem=f"part{part:04}";f=a.native/(stem+".json");values=f.with_suffix(".npz")
    assert sha(f)==h1 and sha(values)==h2
    with np.load(values,allow_pickle=False) as z:pos=z["pos"].tolist()
    assert len(pos)==json.loads(f.read_bytes())["remaining"] and pos==sorted(set(pos))
    if pos:left[str(part)]=pos
count=sum(map(len,left.values()))
assert count==report["totals"]["remaining"]
a.output.write_text(json.dumps({"kind":"unclosed original leaves of the scored threshold lift","target":"0.8372",
    "retained":left,"count":count,"parent_sha256":sha(a.parent),"primary_sha256":sha(a.primary_ledger)},
    separators=(",",":"))+"\n")
print(len(left),count)
