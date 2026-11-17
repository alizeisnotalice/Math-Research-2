import sympy as s,json
from pathlib import Path
L,a,k,d,e,f,u2,u4,eta,rho,G7=s.symbols('L a k d e f u2 u4 eta rho G7')
g2=L**2+u2;g3=-L*g2-k;g4=u4+L**2*g2+2*L*k;g5=eta-e-k*g2-L*g4;g6=rho-f+2*L*e-L*eta+2*L*k*g2+L**2*g4+k**2
P=[1,L,a,L*a+k,d,e,f,1];G=[1,-L,g2,g3,g4,g5,g6,G7]
A7=s.expand(sum(P[j]*G[7-j] for j in range(8)))
U7=a*(eta-e-k*g2)-d*(L*g2+k)
R7=1+L*rho-L**2*eta+2*L**2*k*g2+(L**3+k)*g4+L*k**2+e*(2*L**2+g2)-2*L*f
assert s.expand(A7-G7-U7-R7)==0
S=[0]
for n in range(1,8):S.append(s.expand(-n*P[n]-sum(P[j]*S[n-j] for j in range(1,n))))
want=-7*k*a**2+7*L**2*k*a-7*L*a*d+7*a*e-L**7-7*L**4*k+7*L**3*d-7*L**2*e-7*L*k**2+7*L*f+7*k*d-7
assert s.expand(S[7]-want)==0
out={'status':'passed_formal_A7_and_Newton7','sympy':s.__version__,'U7':str(U7),'R7':str(R7),'s7':str(S[7])};Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
