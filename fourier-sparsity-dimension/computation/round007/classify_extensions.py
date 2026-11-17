#!/usr/bin/env python3
"""Deduplicated span-BFS independent of increasing-basis DFS; classify dimension 15."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'round006'))
from verify_switch_closure import echelon
from analyze_minimum_switches import cover,span


def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    raw=a.input.read_bytes();data=json.loads(raw);K=set(data['K']);began=time.monotonic()
    states={(0,):K};levels=[{'dim':0,'states':1}]
    target=4
    for depth in range(1,target+1):
        nxt={}
        for U,E in states.items():
            unseen=set(E).difference(U)
            while unseen:
                v=min(unseen);coset={v^u for u in U};unseen.difference_update(coset)
                enlarged=tuple(sorted(set(U)|coset))
                if enlarged in nxt:continue
                feasible=E.intersection(v^x for x in E)
                if len(feasible)>=1<<target:nxt[enlarged]=feasible
        states=nxt;levels.append({'dim':depth,'states':len(states)})
    records=[];max_feasible=0
    for U,E in states.items():
        assert len(U)==16 and set(U)<=K
        basis=[x for x,_ in echelon(U).values()]
        G=data['W_generators']+basis;assert len(echelon(G))==15
        minimum=[x for x in span(G) if x.bit_count()==8]
        partition=cover(minimum)
        max_feasible=max(max_feasible,len(E))
        records.append({'quotient_basis':basis,'feasible_extension_count':len(E),'minimum_words':len(minimum),'minimum_span_rank':len(echelon(minimum)),'partition_masks_hex':None if partition is None else [hex(x) for x in partition]})
    assert max_feasible<32 # independent no-five-dimensional-space certificate
    result={'status':'passed','input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'levels':levels,'dimension15_extensions':len(records),'largest_feasible_at_dimension4':max_feasible,'all_have_eight_block_partition':all(r['partition_masks_hex'] is not None for r in records),'records':records,'elapsed_seconds':time.monotonic()-began}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}))

if __name__=='__main__':main()
