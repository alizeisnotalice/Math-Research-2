"""Check the paired Taylor coefficients and their cancellation exactly."""
import json
from pathlib import Path
import sympy as sp

z,c,a,bp,bm,dp,dm=sp.symbols("z c a bp bm dp dm")
out=[]
for n in (4,6,8,10,12,14):
    rows=[]
    for sigma,b,d in ((1,bp,dp),(-1,bm,dm)):
        source=1-sigma*a*z+(a*a-1)*z*z/2
        v=c*z**n*(1+b*z+d*z*z)
        terms=sp.expand(sp.diff(source*v,z,2)/2+sp.diff(source*v*v,z,4)/8)
        rows.append([sp.factor(terms.coeff(z,n-2+j)) for j in range(3)])
    assert rows[0][0]==rows[1][0]
    assert sp.expand((rows[0][1]+rows[1][1])*sp.factorial(n-1)
                     -sp.factorial(n+1)*c*(bp+bm)/2)==0
    assert sp.expand((rows[0][2]-rows[1][2]).subs(bm,-bp)*sp.factorial(n)
                     -sp.factorial(n+2)*c*(dp-dm)/2)==0
    out.append({"n":n,"plus":list(map(str,rows[0])),"minus":list(map(str,rows[1]))})
Path("matched-coefficients.json").write_text(json.dumps({"sympy":sp.__version__,"rows":out},indent=2)+"\n")
