"""Exact tests for the pure sextic escape via one positive diagonal block."""
from fractions import Fraction as F
from math import comb, factorial, lgamma,log,exp
import json,platform
import sympy as sp
N,u,w,b=sp.symbols('N u w b')
t=F(3,8)
lam_sq=(9*t-2)**2/(16*t**3)
threshold_sq=F(4,25)*F(12,5)**3
assert lam_sq==F(121,54)>threshold_sq
p=F(8,19)
assert sum(p**(2*j+1)/F(2*j+1) for j in range(3))>F(4,9)
K=F(18**18*10**10,12**12*3**3*4**4*40**3*10**6)/lam_sq**3
# e<68/25 using exp series through degree 6 and its geometric remainder.
assert sum(F(1,factorial(j)) for j in range(7))+F(1,factorial(7))*F(8,7)<F(68,25)
margin=K/F(68,25)**9
assert margin>1
lam=float(lam_sq)**.5
rows=[]
for N in (100,1000,10000,100000):
    p=3*N//20
    gain=(lgamma(6*p+1)-lgamma(4*p+1)-lgamma(p+1)-p*log(2*N)
        +lgamma(N//2+1)-lgamma(N//2-2*p+1)-2*p*log(lam*N/2))/N
    rows.append({'N':N,'log_block_ratio_over_N':gain})
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,
  'lambda_squared':str(lam_sq),'root_positivity_threshold_squared':str(threshold_sq),
  'exact_stirling_margin_over_rational_e_upper':str(margin),'positive_limiting_excess':(log(float(K))-9)/20,
  'rows':rows},indent=2))
