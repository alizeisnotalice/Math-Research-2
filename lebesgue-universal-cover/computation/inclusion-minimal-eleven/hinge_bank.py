"""Create common integer inner polygons for closed pose and rotation cells."""
from fractions import Fraction as F
import argparse,json,sys
from pathlib import Path
from flint import arb,ctx
from template import build
S=10**7

def q(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

def det(x,y):
    return x[0]*y[1]-x[1]*y[0]

def edges(p):return zip(p,p[1:]+p[:1])

def write_poly(f,p):
    f.write(str(len(p))+' '+' '.join(str(v) for z in p for v in z)+'\n')

def bank(source,count,output,radial=F(".47"),speed=F(".2"),outer=F(".51")):
    ctx.prec=160
    raw=json.loads(source.read_bytes());p=[tuple(F(v,raw['denominator']) for v in z) for z in raw['points']]
    p=sorted(p)
    def half(rows):
        h=[]
        for x,y in rows:
            while len(h)>=2 and (h[-1][0]-h[-2][0])*(y-h[-1][1])<=(h[-1][1]-h[-2][1])*(x-h[-1][0]):h.pop()
            h.append((x,y))
        return h[:-1]
    p=half(p)+half(p[::-1])
    assert len(p)>=3 and all(x*x+y*y<outer**2 for x,y in p)
    assert all(det(a,b)>0 and det(a,b)**2>radial**2*((a[0]-b[0])**2+(a[1]-b[1])**2) for a,b in edges(p))
    # On every differentiable support arc, h'^2+h^2 is at most outer^2.
    # The inradius bound then certifies this angular support speed.
    assert outer**2-radial**2<speed**2
    scale=1-(q(speed)*arb.pi()/count+q(F(2,S)))/q(radial)
    assert scale>0
    p=[(q(x)*S*scale,q(y)*S*scale) for x,y in p]
    def rnd(x):return int((x+q(F(1,2))).floor().unique_fmpz())
    with output.open('w') as f:
        f.write(f'{S} {count}\n')
        for sign in (1,-1):
            for j in range(count):
                s,c=(arb.pi()*(2*j+1)/count).sin_cos()
                v=[(rnd(x*c-sign*y*s),rnd(x*s+sign*y*c)) for x,y in p]
                write_poly(f,v if sign==1 else v[::-1])

def host(c,box,k):
    h,D=build(c,False,k,True,box=box)
    # H includes its approximated sampled disk of radius >.49. Rounding the
    # further contracted vertices by <= sqrt(2)/(2S) is safe since 4*.49>1.
    assert all(det(a,b)>0 and det(a,b)**2*10**4>49**2*D**2*((a[0]-b[0])**2+(a[1]-b[1])**2) for a,b in edges(h))
    v=[tuple(round(F(z*(S-4),D)) for z in t) for t in h]
    sys.path.insert(0,str(c.resolve()));from certify_nine import hull
    return hull(v)

if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='kind',required=True)
    a=sp.add_parser('bank');a.add_argument('source',type=Path);a.add_argument('count',type=int);a.add_argument('output',type=Path)
    a=sp.add_parser('hosts');a.add_argument('comp',type=Path);a.add_argument('boxes',type=Path);a.add_argument('output',type=Path);a.add_argument('--k',type=int,default=192)
    a=p.parse_args();assert __debug__
    if a.kind=='bank':bank(a.source,a.count,a.output)
    else:
        rows=json.loads(a.boxes.read_bytes())
        with a.output.open('w') as f:
            f.write(f'{S} {len(rows)}\n')
            for row in rows:write_poly(f,host(a.comp,row,a.k))
