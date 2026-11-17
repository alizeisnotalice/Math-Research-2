"""Reconstruct every pair, and verify exact positive-weight Cauchy certificate."""
import itertools,json,sys
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).resolve().parent.parent
d=json.loads((p/'n5-ledger.json').read_text());h=d['curves'];idx={s:i for i,s in enumerate(h)};n=d['n'];N=len(h)
assert n==5 and N==32 and len(set(h))==N
assert h==[''.join(x) for x in itertools.product('23',repeat=5)]
def ts(s):
 c=''.join(str(n-int(x)) for x in s);return s,s[::-1],c,c[::-1]
matrix=[[None]*N for _ in h];covered=set();diag=off=0
for j in d['jobs']:
 key=tuple(j['curves']);v=int(j['count']);assert 0<=v
 actual=set(tuple(sorted(x)) for x in zip(ts(key[0]),ts(key[1])))
 assert actual==set(map(tuple,j['members'])) and key==min(actual)
 assert j['diagonal_multiplicity']==sum(a==b for a,b in actual)
 assert j['offdiagonal_multiplicity']==sum(a!=b for a,b in actual)
 diag+=j['diagonal_multiplicity'];off+=j['offdiagonal_multiplicity']
 for a,b in actual:
  assert (a,b) not in covered;covered.add((a,b));i,k=idx[a],idx[b]
  matrix[i][k]=matrix[k][i]=v
assert covered=={(a,b) for i,a in enumerate(h) for b in h[i:]}
assert (diag,off)==(32,496)
C=[matrix[i][i] for i in range(N)];S1=sum(C);S2=sum(map(sum,matrix));A=7063261481362958280000
expected=json.loads((p/'evidence/macbook-aggregate-audit.json').read_text())
assert S1==int(expected['S1']) and S2==int(expected['S2'])
uniform=F(S1*S1,A*S2)
cert=json.loads((p/'weighted-positive-certificate.json').read_text());assert cert['curves']==h
w=list(map(F,cert['weights']));assert all(x>=0 for x in w) and sum(x>0 for x in w)==12
mean=sum(c*x for c,x in zip(C,w))
Mw=[sum(matrix[i][j]*w[j] for j in range(N)) for i in range(N)]
second=sum(x*y for x,y in zip(w,Mw));lower=mean*mean/(A*second)
assert mean==F(cert['mean']) and second==F(cert['second_moment']) and lower==F(cert['lower_bound'])
# Dual inequalities certify optimality within all nonnegative event weights.
assert all(Mw[i]>=C[i] and (not w[i] or Mw[i]==C[i]) for i in range(N))
assert lower>uniform and lower<1
print(json.dumps(dict(n=5,pairs_checked=len(covered),orbit_jobs=len(d['jobs']),S1=str(S1),S2=str(S2),uniform_lower=str(uniform),uniform_lower_decimal=float(uniform),nonnegative_weighted_lower=str(lower),nonnegative_weighted_lower_decimal=float(lower),positive_weights=12,all_32_dual_inequalities_passed=True,scope='Finite n5 rho>=3 lower bounds only'),indent=2))
