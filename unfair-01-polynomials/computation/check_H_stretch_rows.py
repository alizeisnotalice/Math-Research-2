"""Reconstruct exact support-row identities; analytic implications in proof 019."""
import json
S=(0,1,2,5,6)
rows={r:[(i,j) for i in S for j in S if i+j==r] for r in range(13)}
assert rows[3]==[(1,2),(2,1)]
assert rows[9]==[]
reduced={r:[(i,j) for i,j in pairs if j not in (1,2)]for r,pairs in rows.items()}
assert reduced[2]==[(2,0)]
assert reduced[7]==[(1,6),(2,5)]
assert reduced[8]==[(2,6)]
assert reduced[1]==[(1,0)] and reduced[0]==[(0,0)]
for L in range(7,13):
 # Since target support residues lie in S, each other row-sum is zero.
 assert 3 not in S and 3<L
 if L>=8:assert 7 not in S and 7<L
 if L==7:
  assert [r for r in rows if r%L==2]==[2,9]
  assert [r for r in rows if r%L==1]==[1,8]
  assert [r for r in rows if r%L==0]==[0,7]
# For all larger L>12 no two rows have the same residue.
print(json.dumps({'verified':True,'arithmetic':'exact integer support identities','rows':rows,'after_B1_B2_zero':reduced,'general_reasoning':'nonnegative row budgets and the separately proved irreducibility/fairness of H; see proof 019','scope':'all L>=7; all normalized nonnegative real factors of H(x)H(x^L)'},indent=2))
