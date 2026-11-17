#!/usr/bin/env python3
"""Audit C++ full L histograms against independently produced r028 prefix DP."""
from __future__ import annotations
import argparse, json
from fractions import Fraction
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('json',type=Path); ap.add_argument('--reference-n',type=int); ap.add_argument('--reference-dir',type=Path)
    a=ap.parse_args(); d=json.loads(a.json.read_text()); n=d['n']; assert d['u128_checked']
    result={'n':n,'status':'checked','holes':[]}
    reference=None
    if a.reference_n:
        assert a.reference_dir is not None
        reference={x['v']:x for x in json.loads((a.reference_dir/f'tail-n{a.reference_n}.json').read_text())['vertices']}
    for x in d['holes']:
        h=list(map(int,x['histogram_L'])); B=sum(h); A=sum((i+1)*z for i,z in enumerate(h))
        assert B==int(x['B_v_from_D']) and A==int(x['A2_v_from_fibres'])
        item={'v':x['v'],'xy':x['xy'],'B':str(B),'A':str(A),'B_over_A':str(Fraction(B,A))}
        if reference is not None:
            y=reference[x['v']]
            assert h==y['histogram_L'],(x['v'],'histogram')
            assert B==y['B_v_from_D'] and A==y['A2_v_from_fibres'],(x['v'],'totals')
            item['independent_r028_prefix_dp_match']=True
        tails=[sum(h[i:]) for i in range(len(h))]
        ratios=[Fraction(tails[i+1],tails[i]) for i in range(len(tails)-1) if tails[i]]
        item['max_consecutive_tail_ratio']=str(max(ratios,default=Fraction(0)))
        item['tail_ratio_le_4_5_for_this_hole']=all(q<=Fraction(4,5) for q in ratios)
        result['holes'].append(item)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__': main()
