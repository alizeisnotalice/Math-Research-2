"""Exact profiles and root-direction checks for retained ordered ledgers."""
import json,sys
from fractions import Fraction as F
from pathlib import Path

root=Path(__file__).resolve().parent
profiles=[]
for p in [root/'aggregate/raw/n3-aggregate.json',root/'aggregate/raw/n4-aggregate.json',root/'macbook/n5-aggregate-v3.json']:
 d=json.loads(p.read_text());n=d['n'];N=n*n;A=int(d['A2']);D=int(d['D_euclidean_squared_sum']);rows=d['ordered_q_lowb']
 table={(r['q'],r['low_b']):r for r in rows};assert len(table)==len(rows)
 assert sum(int(r['count']) for r in rows)==A
 assert sum(int(r['D_euclidean_squared_sum']) for r in rows)==D
 for r in rows:
  q,b=r['q'],r['low_b'];assert 1<=b<q<N
  s=table[(q,q-b)]
  for field in ('count','D_euclidean_squared_sum'):
   assert b*int(r[field])==(q-b)*int(s[field])
  assert int(r['D_euclidean_squared_sum'])<=2*(n-1)**2*int(r['count'])
 def mass(pred):return sum(int(r['count']) for r in rows if pred(r))
 def moment(pred):return sum(int(r['D_euclidean_squared_sum']) for r in rows if pred(r))
 # Every root direction keeps its exact binomial interleaving weight.
 assert sum((F(r['q'],r['q']-r['low_b'])*int(r['count']) for r in rows),F())==2*A
 eq=F(sum(r['q']*int(r['count']) for r in rows),A)
 eb=F(sum(r['low_b']*int(r['count']) for r in rows),A)
 profiles.append(dict(n=n,nonzero_qb_cells=len(rows),A1=int(d['A1']),A2=A,D2=D,Q=str(F(D,N*A)),R=str(F(A,int(d['A1']))),E_q=str(eq),E_q_over_N=str(eq/N),E_low_b=str(eb),E_low_b_over_N=str(eb/N),p_low_b_1=str(F(mass(lambda r:r['low_b']==1),A)),Q_contribution_low_b_1=str(F(moment(lambda r:r['low_b']==1),N*A)),p_q_le_n=str(F(mass(lambda r:r['q']<=n),A)),root_direction_checks_passed=True))
print(json.dumps(dict(scope='Finite exact original-law profiles only; no asymptotic inference',profiles=profiles),indent=2))
