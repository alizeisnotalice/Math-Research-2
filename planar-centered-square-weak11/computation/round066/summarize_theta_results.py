#!/usr/bin/env python3
"""Exact saved-output aggregation and scope classification; no halo rerun."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def classify(case):
    p=[tuple(map(F,q)) for q in case['points']];m=list(map(F,case['masses']))
    span2=max(max(q[a] for q in p)-min(q[a] for q in p) for a in (0,1))**2
    rows=[]
    for b in case['blocks']:
        score=b['high_unpaid_minus_low_old_credit_normalized'];threshold=b['universal_compensation_threshold']
        passes=F(score['upper'])<=F(threshold['lower'])
        fails=F(score['lower'])>F(threshold['upper'])
        assert passes==b['certifies_universal_compensation_condition']
        fat=all(16*F(v)>9 for v in b['S_width_squared_over_L'])
        comp=next(c for c in b['components'] if c['mask']&b['S_mask'])
        assert comp['mask']&b['S_mask']==b['S_mask']
        rows.append({'mask':b['S_mask'],'rank':b['S_rank'],'L':b['L'],
                     'fat_both_axes':fat,'strong_pass':passes,'strong_strict_fail':fails,
                     'ambiguous':not(passes or fails),'budget_7over2':b['budget_certifies_7over2'],
                     'normalized_score':score,'width_squared_over_L':b['S_width_squared_over_L'],
                     'component_containing_S':comp,'S_is_whole_component':comp['mask']==b['S_mask']})
    components=case['blocks'][0]['components'] if case['blocks'] else []
    return {'id':case['id'],'source_count':len(p),'proper_max_blocks':len(rows),
            'all_proper_maximum_evaluated':case['all_proper_maximum_masks_evaluated'],
            'connected':len(components)==1 if components else None,
            'whole_nonterminal_35_33':33*span2>35*sum(m),
            'all_components_terminal_35_33':all(c['terminal_35_33'] for c in components) if components else None,
            'some_thin_max_block':any(not r['fat_both_axes'] for r in rows),
            'all_max_fat_and_strict_fail':bool(rows) and all(r['fat_both_axes'] and r['strong_strict_fail'] for r in rows),
            'all_max_strict_fail':bool(rows) and all(r['strong_strict_fail'] for r in rows),
            'some_strong_pass':any(r['strong_pass'] for r in rows),'blocks':rows}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--budget-root',type=Path,required=True)
    ap.add_argument('--erosion-root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();assert os.path.ismount('/Volumes/data')
    cases=[];pins=[];ring_stats={'blocks':0,'same_interval_as_actual_high_unpaid':0,'individual_competitors':0,
                               'individual_geometric_ring_upper_at_most_a':0}
    for p in sorted(args.budget_root.glob('stratum-??.json')):
        d=json.loads(p.read_text());assert d['status']=='PASS';cases+=d['cases']
        pins.append({'path':str(p),'sha256':sha(p)})
    ledger=[classify(c) for c in cases];rows=[r for c in ledger for r in c['blocks']]
    for phase in ('pilot','strata'):
        receipt=json.loads((args.erosion_root/f'batch-{phase}.receipt.json').read_text())
        assert receipt['status']=='PASS'
        for job in receipt['jobs']:
            assert job['exit_status']==0 and job['owned_pid_absent']
            p=args.erosion_root/f"{job['job']}-erosion.json";d=json.loads(p.read_text())
            assert sha(p)==job['output_sha256'] and d['status']=='PASS'
            source=args.budget_root/('r64-final-theta.json' if job['job']=='r64-final' else f"{job['job']}.json")
            assert sha(source)==job['budget_result_sha256']==d['budget_result_sha256']
            for c in d['cases']:
                for b in c['blocks']:
                    ring_stats['blocks']+=1
                    ring_stats['same_interval_as_actual_high_unpaid']+=b['ring_union_outside_container']==b['actual_high_unpaid_union_from_budget']
                    ring_stats['individual_competitors']+=len(b['competitors'])
                    ring_stats['individual_geometric_ring_upper_at_most_a']+=sum(r['geometric_ring_upper_at_most_a'] for r in b['competitors'])
    allfail=[c for c in ledger if c['all_max_strict_fail']]
    connected_fail=[c for c in allfail if c['connected'] and c['whole_nonterminal_35_33']]
    fat_fail=[(c['id'],r) for c in ledger for r in c['blocks'] if r['fat_both_axes'] and r['strong_strict_fail']]
    obstruction=[c['id'] for c in connected_fail if c['all_max_fat_and_strict_fail']]
    r64=json.loads((args.budget_root/'r64-final-theta.json').read_text())['cases'][0]['blocks'][0]
    payload={'round':66,'status':'SAVED EXACT OUTPUTS AND RECEIPTS VERIFIED; INDEPENDENT GEOMETRY AUDIT SEPARATE',
             'source_sha256':sha(Path(__file__)),'R63':{'cases':len(cases),'proper_max_blocks':len(rows),
             'budget_pass_blocks':sum(r['budget_7over2'] for r in rows),'strong_pass_blocks':sum(r['strong_pass'] for r in rows),
             'strong_strict_failure_blocks':sum(r['strong_strict_fail'] for r in rows),'ambiguous_blocks':sum(r['ambiguous'] for r in rows),
             'connected_nonterminal_all_max_strict_failure_cases':[c['id'] for c in connected_fail],
             'connected_nonterminal_all_max_fat_and_strict_failure_cases':obstruction,
             'fat_individual_strict_failures':fat_fail,'cases_ledger':ledger},
             'erosion':ring_stats,'R64_final':{k:r64[k] for k in ['S_mask','L','container_area','low_old_credit','high_unpaid_union','low_hole','true_deletion','analytic_envelope_budget','certifies_universal_compensation_condition']},
             'budget_output_pins':pins,
             'scope':'No new geometry by summarizer. Strong failure requires lower(score)>upper(threshold). Thin and disconnected counterexamples are not fat connected-nonterminal obstructions. Equal computed interval bounds are not a separate proof of exact geometric equality.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'blocks':len(rows),'strongpass':sum(r['strong_pass'] for r in rows),'strictfail':sum(r['strong_strict_fail'] for r in rows),'fat_connected_obstruction_count':len(obstruction),'erosion':ring_stats}))


if __name__=='__main__':main()
