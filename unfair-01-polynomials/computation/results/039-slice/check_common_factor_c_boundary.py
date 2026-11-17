#!/usr/bin/env python3
"""Exact boundary example: common Q(x^p) and binary C can coexist.

It validates only the factor/divisibility/binary edges that would remain after
the conjugate-spread lemma.  No F,G/R_F/R_G full slice realization is claimed.
"""
from pathlib import Path
import hashlib,json,sys
import sympy as S
x,y=S.symbols('x y');p,q,L=2,11,6
Q=sum(y**i for i in range(5))
H=S.expand(Q.subs(y,x**p)); U=H
Sslice=S.expand(Q.subs(y,y**L)); V=S.expand(Sslice.subs(y,x**p)); W=S.expand(U*V)
C=S.expand(H.subs(x,x**p)*H.subs(x,x**q))
def binary(P,var=x):
 c=S.Poly(P,var).all_coeffs()
 return set(c)<={0,1}
assert binary(H) and binary(U) and binary(V) and binary(W) and binary(C)
assert S.rem(S.Poly(Sslice,y),S.Poly(Q,y))==0
assert S.rem(S.Poly(V,x),S.Poly(Q.subs(y,x**p),x))==0
assert S.rem(S.Poly(H,x),S.Poly(Q.subs(y,x**p),x))==0
record={'status':'passed_exact_boundary_for_common_factor_plus_C','parameters':{'p':p,'q':q,'L':L,'Q':str(Q)},
 'polynomials':{'H=U':str(H),'S':str(Sslice),'V=S(x^p)':str(V),'W=UV':str(W),'C':str(C)},
 'verified':{'Q_divides_S':True,'Q(x^p)_divides_H':True,'Q(x^p)_divides_U_and_V':True,'H_U_V_W_binary':True,'C_binary':True},
 'scope':'This refutes only an attempted inference from the common-factor conclusion Q(x^p)|H plus binary C/U/V/W to contradiction. It does not supply F,G or q-slice compatibility, so it is not a complete-crossed-model or original-conjecture counterexample.',
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'sympy':S.__version__}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');print(json.dumps(record,indent=2,sort_keys=True))
