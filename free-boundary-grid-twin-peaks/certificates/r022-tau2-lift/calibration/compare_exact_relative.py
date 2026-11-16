#!/usr/bin/env python3
"""Compare r022 strict summaries with existing finite exact certificates."""
import json, math
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 80
base = Path(__file__).resolve().parent
cert = {
  3: base/'inputs/second-birth-survival-n3.json',
  5: base/'inputs/n5-root-moments.json',
  7: base/'inputs/n7-d4-root-moments.json',
}
summ = {
  3: base/'reports/n3-p4096-r64-l16.strict-summary.json',
  5: base/'reports/n5-p8192-r64-l16.strict-summary.json',
  7: base/'reports/n7-p16384-r64-l16.strict-summary.json',
}
out = []
for n in (3,5,7):
    x, s = json.load(open(cert[n])), json.load(open(summ[n]))
    A = int(x['A2'])
    if n == 3:
        D, F, G = int(x['full_euclidean_square_total']), int(x['full_root_F_total']), int(x['full_root_G_total'])
    elif n == 7:
        D, F, G = int(x['squared_euclidean_total']), int(x['firstsum']), int(x['secondsum'])
    else:
        F, G = int(x['firstsum']), int(x['secondsum']); D = (F - 2*G)//4
        assert 4*D == F-2*G
    exact = {'A2_over_factorial': Fraction(A, math.factorial(n*n)), 'Q': Fraction(D, n*n*A), 'V': Fraction(F, 8*A), 'C': Fraction(G, 4*A)}
    row = {'n':n, 'certificate':str(cert[n]), 'A2':str(A),
           'exact':{k:str(Decimal(v.numerator)/Decimal(v.denominator)) for k,v in exact.items()},
           'estimate':{'A2_over_factorial':s['lift_full_normalizer_mean'],'Q':s['Q'],'V':s['V'],'C':s['C']},
           'absolute_error':{k:str(abs(Decimal(s['lift_full_normalizer_mean'] if k=='A2_over_factorial' else s[k])-Decimal(v.numerator)/Decimal(v.denominator))) for k,v in exact.items()}}
    out.append(row)
print(json.dumps(out,indent=2))
