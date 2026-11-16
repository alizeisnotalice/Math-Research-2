#!/usr/bin/env python3
"""Bind separate complete-halo production and compare exact result intervals.

No halo is recomputed here. Owning source and per-result hash are retained.
"""
import hashlib,json,os,sys
from pathlib import Path
from fractions import Fraction as F
out=[]
for name in ['s25-only','r071-positive33','aligned225','contracted-unique']:
 p=Path('computation/results/round072-actual')/(name+'.json')
 if not p.exists():p=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11/computation/results/round072-actual')/(name+'.json')
 d=json.loads(p.read_text());k=['positive_source_count','maximum_mass','maximum_count','chosen_S_is_maximum','full_closure_count','full_eligible_closure_count','eligible_component_count','actual_total_mass','actual_span','nonterminal_gap_span_squared_minus_35total_over33','nonterminal']
 out.append({'name':name,'owner_result_path':str(p),'owner_result_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'full_halo':d['full_halo'],'all_maximum_block_evaluations':[{'chosen_S':a['chosen_S'],'maximum_mask':a['maximum_mask'],'actual_Delta':a['actual_Delta']} for a in d['all_maximum_block_evaluations']]}|{a:d[a] for a in k})
by={r['name']:r for r in out};new=by['contracted-unique'];assert new['maximum_count']==1 and new['eligible_component_count']==1 and new['nonterminal']
nd=new['all_maximum_block_evaluations'][0]['actual_Delta'];od=next(a['actual_Delta'] for a in by['r071-positive33']['all_maximum_block_evaluations'] if a['chosen_S']);assert F(nd['lower'])>F(od['upper']) and F(nd['lower'])>F(3,2) and F(nd['lower'])>F(345909,250000)
result={'python':sys.version,'pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Source-bound reuse of separate full-halo production, not an independent replay; exact Fraction comparisons of intervals only. Independent geometry audit is owned by the geometry review agent.','cases':out,'exact_comparisons':{'contracted_unique_lower_exceeds_r071_chosen_S_upper':True,'contracted_unique_lower_exceeds_3over2':True,'contracted_unique_lower_exceeds_cell_certificate':True,'certified_improvement_over_r071_chosen_S_at_least':str(F(nd['lower'])-F(od['upper']))}}
p=Path('computation/results/round072-cells/complete-geometry-bindings.json');p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'pid':os.getpid(),'exact_comparisons':result['exact_comparisons']},indent=2))
