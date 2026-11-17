#!/usr/bin/env python3
"""Emit the exact 14-worker plan; it starts no process and performs no scan."""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TOTALS = {1: 256, 2: 16384, 3: 262144, 4: 1048576}
WORKERS = {1: 1, 2: 1, 3: 3, 4: 9}  # 14 total workers; caller is the fifteenth process.

def cuts(total: int, pieces: int) -> list[tuple[int,int]]:
    return [(total*i//pieces, total*(i+1)//pieces) for i in range(pieces)]

def main() -> None:
    shards=[]
    for left, workers in WORKERS.items():
        for i,(start,stop) in enumerate(cuts(TOTALS[left],workers)):
            shards.append({"shard_id":f"l{left}-{i:02d}","left_size":left,"start":start,"stop":stop,"edge_subsets":stop-start})
    assert len(shards)==14 and sum(x["edge_subsets"] for x in shards)==sum(TOTALS.values())
    (ROOT/"shard-plan.json").write_text(json.dumps({"N":9,"relation":"R27=h2-5h4+5h6-2h7","workers":14,"families":TOTALS,"shards":shards},indent=2)+"\n")

if __name__ == "__main__": main()
