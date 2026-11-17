"""Analytic Fourier bounds with exact final rational comparisons."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json,platform,hashlib,os,time
import mpmath as mp
start=time.monotonic();mp.mp.dps=80
q=1000
positive=Q(1,2*3**9)*Q(4,11)**1415
unit=Q(3,8)**1500
# Analytic estimates justified in proofs/019-critical-tail-budget.md.
assert positive>8*unit
assert Q(9900,7)<1415
assert 2**700>1000
assert sum(Q(7)**j/factorial(j) for j in range(16))>1000
assert Q(q*q,14)>1500 and Q(q*q,8)>1500
# e<11/4 follows from factorial-series comparison; e>8/3 from positive terms.
assert Q(5,2)+Q(1,4)==Q(11,4)
assert sum(Q(1,factorial(j)) for j in range(4))==Q(8,3)
# AM-GM after allocating x^2/8 to the integrable Gaussian envelope.
assert 4*Q(3,8)*Q(3,2)==Q(3,2)**2
num=lambda z:mp.mpf(z.numerator)/z.denominator
out={'parameters':{'q':q,'x0':3,'core':'0<abs(x-3)<exp(-4): v=1/log(1/abs(x-3))','compact_complement':'abs(x)<=4 outside core: v>=1/4','tails':'abs(x)>=4: v>=3/x^2','target_variance':'M>1'},'exact_comparisons':{'positive_lower':str(positive),'unit_upper':str(unit),'total_negative_upper':str(4*unit),'rotated_real_lower':str(positive/2),'gaussian_modulus_upper':str(unit),'positive_gt_eight_units':True},'display':{'positive_over_unit':str(num(positive)/num(unit)),'log10_rotated_lower':str(mp.log10(num(positive)/2)),'log10_gaussian_upper':str(mp.log10(num(unit)))},'runtime':{'python':platform.python_version(),'mpmath':mp.__version__,'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'elementary analytic kernel estimates plus Fraction final comparisons; not a formalization of all analysis'}
Path(__file__).with_name('round019-critical-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'display':out['display'],'runtime':out['runtime']},indent=2))
