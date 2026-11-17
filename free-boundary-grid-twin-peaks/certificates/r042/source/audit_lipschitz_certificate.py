"""Finite exact 178-curve ledger, nonnegative weights and all dual checks."""
import itertools,json
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).resolve().parent.parent;d=json.loads((p/'n5-lipschitz-ledger.json').read_text());h=d['curves'];N=len(h);idx={s:i for i,s in enumerate(h)}
expected=[''.join(map(str,t)) for t in itertools.product(range(1,5),repeat=5) if all(abs(x-y)<=1 for x,y in zip(t,t[1:]))]
assert h==expected and N==178
def ts(s):
 c=''.join(str(5-int(x)) for x in s);return s,s[::-1],c,c[::-1]
groups={}
for i,a in enumerate(h):
 for b in h[i:]:
  key=min(tuple(sorted(t)) for t in zip(ts(a),ts(b)));z=groups.setdefault(key,[0,0,[]]);z[0 if a==b else 1]+=1;z[2].append((a,b))
assert len(groups)==4134 and len(d['jobs'])==4134
matrix=[[None]*N for _ in h];seen=set();zeros=0
for job in d['jobs']:
 key=tuple(job['curves']);assert key not in seen;seen.add(key);diag,off,members=groups[key]
 assert (job['diagonal_multiplicity'],job['offdiagonal_multiplicity'])==(diag,off)
 v=int(job['count']);assert v>=0
 if job['geometric_zero_certificate']:assert v==0;zeros+=1
 for a,b in members:matrix[idx[a]][idx[b]]=matrix[idx[b]][idx[a]]=v
assert seen==set(groups) and all(v is not None for row in matrix for v in row)
C=[matrix[i][i] for i in range(N)];S1=sum(C);S2=sum(map(sum,matrix));A=7063261481362958280000
old=json.loads((p/'evidence/lipschitz-aggregate-audit.json').read_text());assert S1==int(old['S1']) and S2==int(old['S2'])
cert=json.loads((p/'lipschitz-weighted-positive-certificate.json').read_text());assert cert['curves']==h
w=list(map(F,cert['weights']));assert all(x>=0 for x in w)
Mw=[sum(matrix[i][j]*w[j] for j in range(N)) for i in range(N)]
assert all(Mw[i]>=C[i] and (not w[i] or Mw[i]==C[i]) for i in range(N))
mean=sum(C[i]*w[i] for i in range(N));second=sum(w[i]*Mw[i] for i in range(N));lower=mean*mean/(A*second)
assert mean==F(cert['mean']) and second==F(cert['second_moment']) and lower==F(cert['lower_bound'])
assert lower>F('0.19935098') and lower<=1
print(json.dumps(dict(curves=N,orbit_jobs=len(seen),unordered_pairs=15931,geometric_zero_certificates=zeros,S1=str(S1),S2=str(S2),uniform_lower=str(F(S1*S1,A*S2)),uniform_lower_decimal=float(F(S1*S1,A*S2)),nonnegative_weighted_lower=str(lower),nonnegative_weighted_lower_decimal=float(lower),positive_weights=sum(x>0 for x in w),all_178_dual_inequalities_passed=True,scope='Exact finite n5 rho>=3 event-union lower bounds; no asymptotic inference'),indent=2))
