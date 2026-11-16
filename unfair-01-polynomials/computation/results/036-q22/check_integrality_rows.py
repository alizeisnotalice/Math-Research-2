"""Exact symbolic convolution check for the analytic algebraic-integer proof."""
import hashlib,json,sys
from pathlib import Path
import sympy as S
a,b=S.symbols('a b'); f=S.symbols('f0:44'); g={0:1,2:1,6:a,20:1,22:1,24:1,28:b,42:1}
def row(k,G,F):return S.expand(sum(v*F[k-j] for j,v in G.items() if j<=k))
values={f[0]:1,f[2]:0,f[4]:0,f[6]:1-a,f[8]:a,f[10]:1-a,f[12]:a*a,f[14]:0}
checks={4:0,6:1,8:1,10:1,12:1,14:2*a*a}
actual={k:S.expand(row(k,g,f).subs(values)) for k in checks}
assert all(S.expand(actual[k]-v)==0 for k,v in checks.items())
assert S.diff(row(20,g,f),f[14])==a
assert S.diff(row(22,g,f),f[2])==1
rev={42-j:S.sympify(v).subs(a,1) for j,v in g.items()}
assert S.expand(row(14,rev,f).subs(f[0],1))==f[14]+b
assert S.diff(row(36,rev,f),f[14])==1 and rev[36]==1
out={'status':'passed','claim':'symbolic rows only; analytic integrality argument in proof','python':sys.version,'sympy':S.__version__,'rows':{str(k):str(v) for k,v in actual.items()},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(out,indent=2))
