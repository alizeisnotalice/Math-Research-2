"""Exact resource scan over balanced 5/6 width words, without peak counting."""
from itertools import combinations
from pathlib import Path
import json

from count_protected import HERE, complement, independent_sizes


def resources(h):
    histogram = independent_sizes(''.join(str(int(x)-1) for x in h))
    pairs = list(zip([0]+histogram,histogram))
    layer = max(range(len(pairs)),key=lambda k: sum(pairs[k]))
    return {'profile': h, 'area': sum(map(int,h)), 'protected_vertices': sum(map(int,h))-len(h),
            'independent_sets': sum(histogram), 'largest_layer_states': max(histogram),
            'largest_pair_layer':layer,'dense_payload':8*sum(pairs[layer]),
            'predecessor_lookups':sum(k*s for k,s in enumerate(histogram))}


def main():
    cut_rows=[]
    for ys in combinations(range(11),5):
        xs=['5']*11
        for y in ys: xs[y]='6'
        h=''.join(xs)
        if h>h[::-1]:
            continue
        a,b=resources(h),resources(complement(h))
        key=(max(a['dense_payload'],b['dense_payload']),a['predecessor_lookups']+b['predecessor_lookups'],h)
        cut_rows.append((key,a,b))
    cut_rows.sort()
    expected='55555566666'
    straight=[[resources('55555555555'),resources('66666666666')]]
    out={'scope':'Exact planning only, no regional peak totals',
         'scanned_profile_count_before_reflection':462,
         'scanned_row_reflection_representatives':len(cut_rows),
         'ranking':'first minimize max dense live table bytes; then total predecessor lookups',
         'chosen_pair':[cut_rows[0][1],cut_rows[0][2]],
         'single_turn_pair':[r for r in cut_rows if r[1]['profile']==expected][0][1:],
         'uniform_straight_pair':straight[0],
         'frontier_top_eight':[{'score':list(s),'sides':[a,b]} for s,a,b in cut_rows[:8]]}
    path=HERE/'profile-resource-scan.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print('reflections',len(cut_rows), 'best',cut_rows[0], 'straight', straight[0])


if __name__=='__main__':
    main()
