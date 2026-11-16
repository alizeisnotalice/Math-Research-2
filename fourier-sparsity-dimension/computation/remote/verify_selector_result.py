#!/usr/bin/env python3
"""Independent exact verifier for selector_search JSON outputs."""
import argparse
import hashlib
import json
from pathlib import Path

def rank(mask: int) -> int:
    piv = [0] * 5
    for v in range(32):
        if mask >> v & 1:
            x = v
            for bit in range(4, -1, -1):
                if x >> bit & 1:
                    if piv[bit]:
                        x ^= piv[bit]
                    else:
                        piv[bit] = x
                        break
    return sum(bool(x) for x in piv)

def support(choices: list[int]) -> int:
    if len(choices) != 8 or not 0 <= choices[0] < 4:
        raise ValueError("invalid globally normalized choices")
    bands = [[0] * 8 for _ in range(4)]
    for u, z in enumerate(choices):
        if not 0 <= z < 8:
            raise ValueError("choice out of range")
        bands[z & 3][u] = 1 if u == 0 or z < 4 else -1
    mask = 0
    for c in range(4):
        for a in range(8):
            walsh = sum(bands[c][u] * (-1 if bin(a & u).count("1") & 1 else 1) for u in range(8))
            if walsh:
                mask |= 1 << ((a << 2) | c)
    return mask

parser = argparse.ArgumentParser()
parser.add_argument("result", type=Path)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
raw = args.result.read_bytes()
data = json.loads(raw)
if data["family"] != "signed-selector-u3-y2" or not data["global_sign_normalized"]:
    raise ValueError("wrong family")
start,count=data['start'],data['count']
assert data['total_normalized']==8388608 and 0<=start<=start+count<=8388608
assert sum(r['seen'] for r in data['by_rank'])==count-int(start==0 and count>0)
rows = []
for row in data["by_rank"]:
    if row["seen"]:
        choices = [int(x) for x in row["choices"].split(",")]
        decoded_id=choices[0]*(8**7)+sum(choices[u]*8**(u-1) for u in range(1,8))
        assert decoded_id==row['id'] and start<=decoded_id<start+count
        mask = support(choices)
        # Direct full 32-point truth-table character sums, independent of bands.
        truth=[(1 if choices[u]<4 else -1)*(-1 if ((choices[u]&3)&y).bit_count()%2 else 1)
               for u in range(8) for y in range(4)]
        full=[sum(v*(-1 if (freq&x).bit_count()%2 else 1) for x,v in enumerate(truth))
              for freq in range(32)]
        direct_mask=sum(1 << i for i,v in enumerate(full) if v)
        assert direct_mask==mask
        assert sum(v*v for v in full)==32**2
        if mask != int(row["support_mask_hex"], 16) or bin(mask).count("1") != row["min_s"] or rank(mask) != row["d"]:
            raise AssertionError((row, hex(mask), bin(mask).count("1"), rank(mask)))
        rows.append({"d": row["d"], "s": row["min_s"], "ratio": row["min_s"] / row["d"] ** 2, "id": row["id"]})
args.out.write_text(json.dumps({"status": "passed", "result_sha256": hashlib.sha256(raw).hexdigest(),
                              "coverage_counts_checked":True,"witness_ids_checked":True,
                              "full_truth_table_walsh_checked":True,
                              "scope":"checks witnesses and counts; does not independently re-enumerate the family",
                              "verified_minima": rows}, indent=2) + "\n")
