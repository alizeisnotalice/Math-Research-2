"""Symbolic coefficient checks for the S234 odd/even low-row reduction."""
import sympy as s,json
from pathlib import Path
L,a,k,d,e,f,u2,u4,eta,G6=s.symbols('L a k d e f u2 u4 eta G6')
P=[1,L,a,L*a+k,d,e,f,1];g2=L**2+u2;g3=-L*g2-k;g4=u4+L**2*g2+2*L*k;g5=eta-e-k*g2-L*g4;G=[1,-L,g2,g3,g4,g5,G6]
A4=s.expand(sum(P[j]*G[4-j] for j in range(5)))
assert s.expand(A4-d-a*u2-u4)==0
A5=s.expand(sum(P[j]*G[5-j] for j in range(6)))
A6=s.expand(sum(P[j]*G[6-j] for j in range(7)))
U6=a*(u4+L*k)+d*(u2+L**2)
R6=G6+f-2*L*e+L*eta-2*L*k*g2-L**2*g4-k**2
assert s.expand(A5-(eta-k*a-L*d))==0 and s.expand(A6-U6-R6)==0
S=[0]
for n in range(1,7):S.append(s.expand(-n*P[n]-sum(P[j]*S[n-j] for j in range(1,n))))
assert s.expand(S[5]-(5*k*a+5*L*d-L**5-5*L**2*k-5*e))==0
want=-2*a**3-6*L*k*a+6*a*d-6*L**2*d+L**6+6*L**3*k+3*k**2+6*L*e-6*f
assert s.expand(S[6]-want)==0
out={'status':'passed_formal_A4_A5_A6_and_Newton56','sympy':s.__version__,'A5':str(A5),'A6':str(A6),'rho6':str(R6),'U6':str(U6),'s5':str(S[5]),'s6':str(S[6])};Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
