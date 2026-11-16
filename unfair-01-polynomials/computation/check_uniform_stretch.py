"""Exact row checks for uniform Hd stretch theorem; no CAS or float use."""
import json

def check(d):
 S={0,1,2,d-1,d}
 def row(n):return sorted((i,n-i)for i in S if n-i in S)
 assert row(3)==[(1,2),(2,1)]
 assert row(d+1)==[(1,d),(2,d-1),(d-1,2),(d,1)]
 assert row(d+3)==[]
 def reduced(n):return [(i,j)for i,j in row(n)if j not in (1,2)]
 assert reduced(d+2)==[(2,d)]
 assert reduced(d+1)==[(1,d),(2,d-1)]
 assert reduced(2)==[(2,0)]
 L=d+1
 assert 2+2*L>2*d
 assert [r for r in range(2*d+1)if r%L==2]==[2,d+3]
 # Both stages of the proof use residues outside S, not arbitrary coefficients.
 assert 3 not in S and 3<L
 return True
if __name__=='__main__':
 assert all(check(d)for d in range(6,501))
 assert 8 in {i+j for i in (0,1,2,4,5)for j in (0,1,2,4,5)}
 print(json.dumps({'verified_d_range':[6,500],'exact_row_checks':True,'d5_boundary_obstruction':8,'general_proof':'explicit sumset classification plus nonnegative row arguments; finite instances are supplementary'},indent=2))
