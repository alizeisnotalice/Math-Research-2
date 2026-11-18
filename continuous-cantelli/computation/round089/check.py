import sympy as s
v,c=s.symbols('v c')
w=1-v+c*v**2
level=1/v-2*s.log(w)
print('edge',s.series(level,v,0,3))
print('latent',s.series(w*level/(level-1),v,0,3))
assert s.series(w*level/(level-1),v,0,3).removeO().expand()==1+c*v**2
