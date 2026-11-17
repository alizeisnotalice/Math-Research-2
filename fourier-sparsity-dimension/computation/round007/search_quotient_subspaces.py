#!/usr/bin/env python3
"""Complete canonical DFS for a target-dimensional subspace inside a finite K."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def find_space(K,target):
    threshold=1<<target
    stats={'nodes':0,'extensions_considered':0,'cardinality_prunes':0,'largest_depth':0}
    def dfs(span,feasible,basis,lower):
        stats['nodes']+=1;stats['largest_depth']=max(stats['largest_depth'],len(basis))
        if len(basis)==target:return basis
        if len(feasible)<threshold:
            stats['cardinality_prunes']+=1;return None
        # The increasing minimum-outside-span basis of any final subspace is unique.
        # Each selected vector must be the least member of its current span coset.
        candidates=[v for v in sorted(feasible) if v>lower and v and min(v^s for s in span)==v]
        for v in candidates:
            stats['extensions_considered']+=1
            nxt=feasible.intersection(x^v for x in feasible)
            if len(nxt)<threshold:
                stats['cardinality_prunes']+=1;continue
            got=dfs(span+[v^s for s in span],nxt,basis+[v],v)
            if got is not None:return got
        return None
    return dfs([0],set(K),[],0),stats


def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--target',type=int,default=5);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    raw=a.input.read_bytes();d=json.loads(raw);K=set(d['K']);assert 0 in K
    start=time.monotonic();basis,stats=find_space(K,a.target)
    result={'status':'completed','target_dimension':a.target,'target_exists':basis is not None,'witness_basis':basis,'K_size':len(K),'stats':stats,'elapsed_seconds':time.monotonic()-start,'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'exact finite K only'}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
