"""Check the fourth-order saddle-centered raw/delivered cancellation."""
import json
from pathlib import Path
import platform
import sympy as s

h,w=s.symbols("h w")
def jet(f):
    return s.series(f,h,0,5).removeO().expand()
Q=1+h**2/2
a=w*h-w*h**3
delta=jet(-Q*s.log(1+a))
sqrt=jet(s.sqrt(1+h*h*delta))
advance=jet(w*h+h*h*delta+w*h*(sqrt-1))
assert jet(advance-a)==0
den=jet(1+Q*h*h/(1+a)*(1+w*h/(2*sqrt)))
amp=jet(jet(s.exp(delta))*sqrt*jet((1+Q*h*h)/den))
expected=1-w*h+w*w*h*h+(w/2-w**3)*h**3+(w**4-s.Rational(3,2)*w*w)*h**4
assert s.simplify(amp-expected)==0
assert s.Rational(3,4)-s.Rational(3,2)*s.Rational(1,2)==0
out=dict(implicit_a=str(a),delta=str(delta),prefactor=str(amp),second_gaussian_moment="1/2",fourth_gaussian_moment="3/4",
         fourth_order_average="0",python=platform.python_version(),sympy=s.__version__)
Path('cancellation.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
