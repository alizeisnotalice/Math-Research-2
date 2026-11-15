#!/usr/bin/env python3
"""Metadata/exact ledger summary from frozen actual outputs, no halo rerun."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    assert os.path.ismount('/Volumes/data')
    data=json.loads((args.root/'structural-inputs.json').read_text())
    inp=[r for s in data['strata'] for r in s['inputs']]
    def identity(r):return tuple(sorted((x,y,m) for (x,y),m in zip(r['points_quarters'],r['mass_sixteenths'])))
    results={};pins=[]
    for kind in ('deletion-batch','allcuts-batch'):
        receipt=json.loads((args.root/kind/'batch-receipt.json').read_text())
        rows=[]
        for job in receipt['jobs']:
            path=args.root/kind/f"stratum-{job['stratum']:02}.json"
            d=json.loads(path.read_text());assert d['status']=='PASS'
            assert sha(path)==job['output_sha256'] and job['exit_status']==0 and job['owned_pid_absent']
            assert d['script_sha256']==job['source_sha256']=='5e90e5f06926b6716f56c2194b0cedeadf15d761b85173e3829d707ccc8cdf36'
            assert d['input_sha256']==sha(args.root/'structural-inputs.json')
            pins.append({'source':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
            rows+=d['cases']
        assert sorted(r['input']['id'] for r in rows)==sorted(r['id'] for r in inp)
        results[kind]={'receipt':receipt,'rows':rows}
    deletion=results['deletion-batch']['rows'];cuts=results['allcuts-batch']['rows']
    maxleaders=sorted([r for r in deletion if r['deletions']['max_mass_min_upper'] is not None],key=lambda r:F(r['deletions']['max_mass_min_upper']),reverse=True)
    leaders=[]
    for r in maxleaders[:3]:
        d=r['deletions'];p=r['input']['points_quarters'];m=r['input']['mass_sixteenths'];D2=1<<(2*r['bits'])
        blocks=[b for b in d['blocks'] if 'max-mass-strictly-eligible' in b['selection_reasons']]
        winner=min(blocks,key=lambda b:F(b['ratio_upper']));mask=winner['mask'];ids=[i for i in range(len(m)) if mask>>i&1]
        widths=[max(p[i][a] for i in ids)-min(p[i][a] for i in ids) for a in (0,1)];L=sum(m[i] for i in ids)
        leaders.append({'id':r['input']['id'],'terminal':r['all_components_terminal_35_33'],
                        'mask':mask,'rank':len(ids),'mass':str(F(L,16)),
                        'widths':list(map(lambda w:str(F(w,4)),widths)),
                        'width_squared_over_L':[str(F(w*w,L)) for w in widths],
                        'both_axes_hard_3over4':all(16*w*w>9*L for w in widths),
                        'minimum_maxmass_upper':winner['ratio_upper'],
                        'minimum_maxmass_lower':winner['ratio_lower'],
                        'full_outer_area':str(F(int(r['full_outer_area_numerator']),D2)),
                        'remaining_halo_lower_area':str(F(int(winner['remaining_inner_area_numerator']),D2)),
                        'old_halo_overlap_credit_lower':str(F(int(r['full_outer_area_numerator'])-int(winner['loss_numerator_upper']),D2))})
    summary={'round':63,'status':'ALL SAVED RECEIPTS AND SOURCE IDENTITIES VERIFIED; INDEPENDENT LEADER REPLAY SEPARATE',
             'source_sha256':sha(Path(__file__)),'input_sha256':sha(args.root/'structural-inputs.json'),
             'cases_per_task':len(inp),'order_free_distinct_inputs':len({identity(r) for r in inp}),
             'source_counts':[16,17,18,19,20,20],'strata':12,
             'deletion':{'actual_deletion_rank':1,'rank3_evaluated':False,
                         'cases':len(deletion),'terminal_cases':sum(r['all_components_terminal_35_33'] for r in deletion),
                         'nonterminal_cases':sum(not r['all_components_terminal_35_33'] for r in deletion),
                         'all_maximum_mass_ties':sum(len(r['deletions']['all_maximum_strictly_eligible_masks']) for r in deletion),
                         'improper_full_input_ties':sum(r['deletions']['maximum_strictly_eligible_full_input'] for r in deletion),
                         'proper_maximum_ties_evaluated':sum('max-mass-strictly-eligible' in b['selection_reasons'] for r in deletion for b in r['deletions']['blocks']),
                         'deduplicated_block_evaluations':sum(len(r['deletions']['blocks']) for r in deletion),
                         'maximum_block_leaders':leaders,'receipt_wall_seconds':results['deletion-batch']['receipt']['wall_seconds']},
             'allcuts':{'cases':len(cuts),'cuts':sum(r['all_cuts'] for r in cuts),
                        'compact_cuts':sum(r['compact_cuts'] for r in cuts),'noncompact_cuts':sum(r['noncompact_cuts'] for r in cuts),
                        'terminal_improved_cuts':sum(r['terminal_improved_cuts'] for r in cuts),
                        'all_certify_B':all(r['certifies_B'] for r in cuts),
                        'max_geometric_bound':str(max(F(r['geometric_bound']) for r in cuts)),
                        'receipt_wall_seconds':results['allcuts-batch']['receipt']['wall_seconds']},
             'source_output_pins':pins,'scope':'Deletion and allcut tasks separate. No rank3 execution. No new general C2 bound; no repeat72 replay.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('cases_per_task','order_free_distinct_inputs')}))


if __name__=='__main__':main()
