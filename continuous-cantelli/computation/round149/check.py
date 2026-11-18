import sympy as s,json
x,z=s.symbols('x z')
E=s.exp(z*x-z*z/2)
row=-z*z*E/2
for j in range(2,11):
    r=s.diff(row,z,j).subs(z,0).expand()
    assert s.Poly(r,x).LC()==-s.binomial(j,2)
for w in [1,1+s.I]:
    for j in range(7):
        r=s.expand(s.diff(row,z,j)/E).subs(z,w).expand()
        assert s.simplify(s.Poly(r,x).LC()+w**2/2)==0
print(json.dumps(dict(sympy=s.__version__,zero_orders=list(range(2,11)),other_orders=list(range(7)),parameters=['1','1+i'])))
