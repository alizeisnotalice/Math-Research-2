"""Symbolic coefficient audit for the two-sided q-slice identities."""
from pathlib import Path
import hashlib,json,os,sys
import sympy as S

assert os.path.ismount('/Volumes/data')
x,y=S.symbols('x y')
p,q=2,5
F=1+3*x-2*x**2+5*x**5-7*x**7+x**10
G=1-4*x+6*x**3+2*x**5-x**6+3*x**10
def sliceq(P):
    return sum(P.coeff(x,q*j)*y**j for j in range(S.degree(P,x)//q+1))
RF,RG,H0=sliceq(F),sliceq(G),sliceq(S.expand(F*G))
A=S.expand(F.subs(x,x**p)*G.subs(x,x**q))
B=S.expand(F.subs(x,x**q)*G.subs(x,x**p))
C=S.expand((F*G).subs(x,x**p)*(F*G).subs(x,x**q))
U=S.expand(F*RG.subs(y,x**p));V=S.expand(G*RF.subs(y,x**p))
T=S.expand((F*G)*H0.subs(y,x**p))
assert all(A.coeff(x,q*k)==V.coeff(x,k) for k in range(30))
assert all(B.coeff(x,q*k)==U.coeff(x,k) for k in range(30))
assert all(C.coeff(x,q*k)==T.coeff(x,k) for k in range(50))
# Full residue slices underlying the stated missing step.
def residue(P,r):
    return sum(P.coeff(x,r+q*j)*y**j for j in range((S.degree(P,x)-r)//q+1) if r+q*j<=S.degree(P,x))
for r in range(q):
    Fr,Gr,Hr=residue(F,r),residue(G,r),residue(S.expand(F*G),r)
    Apiece=S.expand(G.subs(x,y)*Fr.subs(y,y**p))
    Bpiece=S.expand(F.subs(x,y)*Gr.subs(y,y**p))
    Cpiece=S.expand((F*G).subs(x,y)*Hr.subs(y,y**p))
    assert all(A.coeff(x,p*r+q*k)==Apiece.coeff(y,k) for k in range(40))
    assert all(B.coeff(x,p*r+q*k)==Bpiece.coeff(y,k) for k in range(40))
    assert all(C.coeff(x,p*r+q*k)==Cpiece.coeff(y,k) for k in range(70))
assert S.expand(H0-(RF*RG+y*sum(residue(F,r)*residue(G,q-r) for r in range(1,q))))==0
# Cross-residue terms need not vanish formally.
F2,G2=1+x,1+x**2
def slice3(P):
    return sum(P.coeff(x,3*j)*y**j for j in range(S.degree(P,x)//3+1))
assert slice3(F2)*slice3(G2)==1 and slice3(S.expand(F2*G2))==1+y

def slice_n(P,n):
    return sum(P.coeff(x,n*j)*y**j for j in range(S.degree(P,x)//n+1))
def binary(P):
    return all(c in (0,1) for c in S.Poly(S.expand(P),x).all_coeffs())

# A complete long crossed model where H_0 != R_F R_G; all claims use exact Z[x].
pL,qL=2,3
FL=S.expand(sum((-1)**k*x**(6*k) for k in range(19))*(1+x**16))
GL=S.expand((1+x**4)*(1+x**6)*(1+x**9)*(1+x**446))
HL=S.expand((1+x**114)*(1+x**4)*(1+x**9)*(1+x**16)*(1+x**446))
AL=S.expand(FL.subs(x,x**pL)*GL.subs(x,x**qL))
BL=S.expand(FL.subs(x,x**qL)*GL.subs(x,x**pL))
CL=S.expand(HL.subs(x,x**pL)*HL.subs(x,x**qL))
AT=S.expand((1+x**228)*(1+x**18)*(1+x**27)*(1+x**32)*(1+x**1338))
BT=S.expand((1+x**342)*(1+x**8)*(1+x**12)*(1+x**48)*(1+x**892))
RFL,RGL,H0L=slice_n(FL,qL),slice_n(GL,qL),slice_n(HL,qL)
assert S.expand(FL*GL-HL)==0 and S.expand(AL-AT)==0 and S.expand(BL-BT)==0
assert binary(HL) and binary(AL) and binary(BL) and binary(CL)
assert len(S.Poly(CL,x).terms())==1024 and S.degree(CL,x)==2945
assert S.expand(H0L-RFL*RGL)!=0
assert all(c>=0 for c in S.Poly(S.expand(H0L-RFL*RGL),y).all_coeffs())

# Exact audit of the universal signed-mother family at p=3,q=5,b=11.
pU,qU,bU=3,5,11
FU=S.expand(sum((-1)**k*x**(pU*qU*k) for k in range(bU)))
GU=S.expand((1+x**(pU**2))*(1+x**(pU*qU))*(1+x**(qU**2)))
HU=S.expand((1+x**(pU**2))*(1+x**(qU**2))*(1+x**(pU*qU*bU)))
AU=S.expand(FU.subs(x,x**pU)*GU.subs(x,x**qU))
BU=S.expand(FU.subs(x,x**qU)*GU.subs(x,x**pU))
CU=S.expand(HU.subs(x,x**pU)*HU.subs(x,x**qU))
AUT=S.expand((1+x**(pU**2*qU*bU))*(1+x**(pU*qU**2))*(1+x**(qU**3)))
BUT=S.expand((1+x**(pU*qU**2*bU))*(1+x**(pU**3))*(1+x**(pU**2*qU)))
assert S.expand(FU*GU-HU)==0 and S.expand(AU-AUT)==0 and S.expand(BU-BUT)==0
assert binary(HU) and binary(AU) and binary(BU) and binary(CU)
assert -1 in S.Poly(FU,x).all_coeffs()
assert S.degree(HU,x)==199 and len(S.Poly(HU,x).terms())==8
assert S.degree(CU,x)==1592 and len(S.Poly(CU,x).terms())==64

# Exact recurrence audit for the integrality-transfer lemma, including gcd(p,q)>1.
pI,qI=4,6
FI=1-2*x+3*x**2-x**5+2*x**7
GI=1+4*x-3*x**3+x**4
HI=S.expand(FI*GI)
AI=S.expand(FI.subs(x,x**pI)*GI.subs(x,x**qI))
N=max(S.degree(FI,x),S.degree(GI,x),S.degree(HI,x),S.degree(AI,x))
frecovered=[S.Integer(1)]
for k in range(1,S.degree(FI,x)+1):
    d=S.gcd(pI,qI)
    left_old=sum(frecovered[k-(qI//d)*ell]*HI.coeff(x,(pI//d)*ell)
                 for ell in range(1,k//(qI//d)+1))
    right=sum(AI.coeff(x,pI*k-qI*j)*frecovered[j]
                  for j in range((pI*k)//qI+1))
    frecovered.append(S.expand(right-left_old))
assert all(frecovered[k]==FI.coeff(x,k) for k in range(len(frecovered)))
grecovered=[S.Integer(1)]
for k in range(1,S.degree(GI,x)+1):
    grecovered.append(S.expand(HI.coeff(x,k)-sum(FI.coeff(x,i)*grecovered[k-i] for i in range(1,k+1))))
assert all(grecovered[k]==GI.coeff(x,k) for k in range(len(grecovered)))

out={'status':'passed','p':p,'q':q,'sympy':S.__version__,'A_projection':str(V),'B_projection':str(U),'C_projection':str(T),'cross_residue_example':{'RF_RG':'1','H0':'1+y'},'long_crossed_H0_not_product':{'p':pL,'q':qL,'C_terms':1024,'C_degree':2945},'universal_signed_family_instance':{'p':pU,'q':qU,'b':bU,'H_terms':8,'C_terms':64,'C_degree':1592},'integrality_transfer_noncoprime_regression':{'p':pI,'q':qI,'F_degree':int(S.degree(FI,x)),'G_degree':int(S.degree(GI,x))},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version}
Path(__file__).with_name('projection-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
