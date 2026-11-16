#!/usr/bin/env python3
"""Exact arithmetic crosschecks of round 9 archived ledgers."""
from pathlib import Path
from fractions import Fraction
from math import factorial
import json
R=Path(__file__).resolve().parent
C=R/'certificates/r009-n9'
def read(p):return json.loads(p.read_text())
a=read(C/'n9-result.json');h=int(a['h']);H=int(a['p_numerator']);T=int(a['T'])
assert T==64*63*H
assert 16*H-4*h==1845019661294621338455581712717756278071765125592476434391513773922680832000>0
cpp=read(C/'n9-reflection-layer-audit-cpp.json');py=read(C/'n9-reflection-layer-audit-python.json')
for row in cpp['layers']:
 k=row['k']
 for key in ['raw','fixed','reps']:assert row[key]==py[key][k]
assert cpp['raw_total']==py['raw_total']==877160008
assert cpp['rep_total']==py['rep_total']==438672693
assert cpp['fixed_total']==py['fixed_total']==185378
assert 48*factorial(64)<2**302
assert 877160008*14*factorial(64)<2**330
n8=read(C/'n8-corrected-crosscheck.json');old=read(R/'results/n8-kernel-counterexample-python-r008.json')['n8']
assert (int(n8['h']),int(n8['p_numerator']),int(n8['T']))==(old['hS'],old['H'],old['pair_total'])
for n in [2,3]:
 d=read(R/f'results/minimum-deletion-n{n}-r009.json');assert d['A1_direct']==d['A1_punctured_sum']
 assert d['B1_punctured_nonadjacent_sum']+d['root_adjacent_minimum_count']==d['A1_direct']
 d=read(R/f'results/original-law-tree-n{n}-r009.json')
 assert d['proposal_terminal_mass']=='1'
 assert sum(x['A2'] for x in d['forced_second_birth_strata'].values())==d['A2_weighted']
 assert sum(x['euclidean_squared_total'] for x in d['forced_second_birth_strata'].values())==d['euclidean_squared_total_weighted']
assert read(R/'results/inert-fibres-r009.json')['maximum_fibre']==4
print(json.dumps({'status':'passed','n9_scaled':str(Fraction(16*H,h)),'n9_scaled_decimal':float(Fraction(16*H,h)),'layer_comparison':'all raw/fixed/reps agree','original_law_full_trees':[2,3]},indent=2))
