#!/usr/bin/env python3
"""Read-only consistency diagnostic for two finite n=24 SMC summaries.

It does not estimate sampling uncertainty and does not pool proposal methods.
"""
from __future__ import annotations
import json
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec=80

src=Path(__file__).resolve().parent.parent / 'r025-genealogy' / 'n24-genealogy-diagnostic.json'
x=json.loads(src.read_text())
out={'scope':'Plug-in checks of necessary exact identities using finite, non-pooled SMC ratio summaries; no CI, significance, or implementation-failure claim.' ,'runs':{}}
for label,run in x['runs'].items():
    z=run['terminal_weight_diagnostic']
    mb=Decimal(z['E_tau2_L'])
    hist=[Decimal(v) for v in z['complete_K2_L_weighted_distance_histogram']]
    assert sum(hist) > Decimal('0.999999') and sum(hist) < Decimal('1.000001')
    b_over_a=Decimal(1)/mb
    linear=sum((Decimal(d)/(2*(d-1))*hist[d] for d in range(2,len(hist))),Decimal())
    coeffsum=sum((Decimal(d-2)/(2*(d-1))*hist[d] for d in range(2,len(hist))),Decimal())
    row={'mbar_plugin':str(mb),'B_over_A_plugin':str(b_over_a),
         'histogram_total':str(sum(hist)),'linear_rhs_plugin':str(linear),
         'linear_slack_BA_le_rhs':str(linear-b_over_a),
         'equivalent_slack_1_BA_minus_coeffsum':str(Decimal(1)-b_over_a-coeffsum),
         'thresholds':[]}
    for r in (2,6,12,18):
        tail=sum(hist[r+1:])
        upper=Decimal(2*r)/(r-1)*(Decimal(1)-b_over_a)
        row['thresholds'].append({'r':r,'tail_plugin':str(tail),'necessary_upper_plugin':str(upper),
                                  'slack_upper_minus_tail':str(upper-tail)})
    out['runs'][label]=row
Path(__file__).with_name('n24-finite-consistency.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
