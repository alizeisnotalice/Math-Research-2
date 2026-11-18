#!/usr/bin/env python3
"""Re-create the exact Metal input from the archived four terminal profiles."""
import json
import struct
from pathlib import Path

rows=json.loads(Path("candidates.json").read_text())
n=max(len(r["offset"]) for r in rows)
with Path("profiles.bin").open("wb") as f:
    f.write(struct.pack("<II",len(rows),n))
    for row in rows:
        v=[row["R"],row["dx"],len(row["offset"]),row["case"],row["layer"]]+row["offset"]
        v+=[0.]*(n+5-len(v))
        f.write(struct.pack("<"+"f"*len(v),*v))
