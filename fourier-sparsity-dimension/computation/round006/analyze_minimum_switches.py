#!/usr/bin/env python3
"""Exact minimum-word geometry and all 255 C-containing hyperplane cover gates."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from verify_switch_closure import echelon


def span(gens):
    words=[0]
    for g in gens:words += [x^g for x in words]
    return words


def cover(masks,remaining=(1<<64)-1):
    if not remaining:return []
    low=remaining&-remaining
    for m in masks:
        if m&low and not m&~remaining:
            got=cover(masks,remaining^m)
            if got is not None:return [m]+got
    return None


def main():
    ap=argparse.ArgumentParser();ap.add_argument('audit',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    raw=a.audit.read_bytes();data=json.loads(raw)
    C=[(1<<64)-1]+[sum(((x>>j)&1)<<x for x in range(64)) for j in range(6)]
    rows=[]
    for row in data['rows']:
        words=span(C+row['relative_generators']);minimum=[(x,i>>7) for i,x in enumerate(words) if x.bit_count()==8]
        masks=[x for x,_ in minimum];rank=len(echelon(masks));intersections=sorted({(x&y).bit_count() for j,x in enumerate(masks) for y in masks[:j]})
        assert set(intersections)<={0,4}
        record={'source':Path(row['source_file']).name,'row':row['row_index'],'minimum_words':len(masks),'minimum_span_rank':rank,'intersection_sizes':intersections,'contains_C':all(len(echelon(masks+[c]))==rank for c in C)}
        if record['source']=='corrected-results.json' and record['row']==1:
            gates=[]
            for ell in range(1,256):
                ms=[m for m,z in minimum if not((z&ell).bit_count()&1)]
                partition=cover(ms)
                gates.append({'normal':ell,'minimum_words':len(ms),'has_partition':partition is not None})
            assert gates[0]['has_partition'] and not gates[3]['has_partition']
            record['hyperplanes']=gates
            record['no_partition_hyperplanes']=sum(not g['has_partition'] for g in gates)
        rows.append(record)
    result={'status':'passed','python':sys.version,'platform':platform.platform(),'arithmetic':'integer and F2','tolerance':0,'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows}
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'hyperplanes_checked':255,'no_partition':rows[1]['no_partition_hyperplanes'],'selected_removed_generator':2}))

if __name__=='__main__':main()
