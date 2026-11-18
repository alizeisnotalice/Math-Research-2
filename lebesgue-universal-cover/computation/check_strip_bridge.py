"""Check strip-code examples and the area budget, using directed Arb intervals."""
from fractions import Fraction
import json
from flint import arb, ctx
ctx.prec = 160

def num(q):
    q = Fraction(q)
    return arb(q.numerator)/q.denominator

def scale(m):
    return (1+num(Fraction(1,m*m)))/(arb.pi()/(2*m)).cos()

def check(m,n,phase):
    angle = arb.pi()/m
    amp = num(Fraction(1,2*(n*n-1)))
    # g(0)=g'(0)=0 after translating the initial antipodal midpoint to zero.
    def g(t):
        return amp*((n*t+phase).cos()-phase.cos()*t.cos()+n*phase.sin()*t.sin())
    k=[]
    for j in range(m):
        v=m*m*g(angle*j)
        nearest=round(float(v.mid()))
        assert abs(v-nearest)<num(Fraction(1,2))
        k.append(nearest)
    assert k[0]==0 and abs(k[1])<=2
    for a,b,c in zip(k,k[1:],k[2:]):
        assert abs(a+c-2*angle.cos()*b)<7
    return len(set(k))

examples=[check(m,n,num(Fraction(3,17))) for m in (5,16,32,64) for n in (3,5,9)]
assert all(scale(m)**2-1 < num(Fraction(6,m*m)) for m in (4,8,16,32,64,128,256))
print(json.dumps(dict(precision=ctx.prec,normalized_example_count=len(examples),
    relative_area_error={str(m):str(scale(m)**2-1) for m in (16,32,64,128,256)}),indent=2))
