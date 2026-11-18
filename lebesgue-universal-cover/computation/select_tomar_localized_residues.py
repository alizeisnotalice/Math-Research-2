"""Build an exact next-stage original-leaf list from every certified higher-localization complement receipt."""
import argparse
import hashlib
import json
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
for key in ("previous_selection","receipts","next_selection"):
    p.add_argument(key,type=Path)
a=p.parse_args()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
previous=json.loads(a.previous_selection.read_bytes());expected=previous["retained"]
frozen={};parents={}
assert {f.stem for f in a.receipts.glob("part*.json")}=={f"part{int(k):04}"for k in expected}
for k,original in sorted(expected.items(),key=lambda kv:int(kv[0])):
    f=a.receipts/f"part{int(k):04}.json";r=json.loads(f.read_bytes())
    assert r["selection_sha256"]==sha(a.previous_selection) and r["part"]==int(k)
    assert r["tested"]==len(original) and r["closed"]+len(r["remaining"])==r["tested"]
    rem=r["remaining"];assert rem==sorted(set(rem)) and set(rem)<=set(original)
    if rem:frozen[k]=rem
    parents[k]=sha(f)
count=sum(map(len,frozen.values()))
result={"kind":"outstanding original leaves after lower-target area checks",
        "target":previous["target"],"totals":[count,0],"retained":frozen,
        "prior_selection_sha256":sha(a.previous_selection),"prior_receipt_sha256":parents}
a.next_selection.write_text(json.dumps(result,separators=(",",":"))+"\n")
print(len(frozen),count)
