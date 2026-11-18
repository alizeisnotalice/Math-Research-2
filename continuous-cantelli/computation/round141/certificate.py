"""Certified central first-resolvent derivative via dyadic boxes, Python stdlib."""
import json
import platform
from decimal import Decimal, getcontext
from functools import cache

getcontext().prec = 96
S = 2**64
DX = 512
DZ = 512
XMAX = 12
ZMAX_N = 9*DZ//16
def ceildiv(a,b): return -((-a)//b)
def down(a,b): return a*b//S
def up(a,b): return ceildiv(a*b,S)
def divlo(a,b): return a*S//b
def divhi(a,b): return ceildiv(a*S,b)
def ratlo(a,b): return a*S//b
def rathi(a,b): return ceildiv(a*S,b)
def encode(value, upper):
    n,d=value.as_integer_ratio()
    return ceildiv(n*S,d) if upper else n*S//d
@cache
def trans(op,n,upper):
    # Integer n/S has at most 64 fractional digits and fewer than 10
    # integer digits at all calls, hence Decimal construction is exact.
    arg=Decimal(n)/Decimal(S)
    val=arg.exp() if op=="exp" else arg.ln()
    enclosure=val.next_plus() if upper else val.next_minus()
    return encode(enclosure,upper)
def elo(n): return max(0,trans("exp",n,False))
def ehi(n): return trans("exp",n,True)
# Archimedes' rational bounds 333/106 < pi < 355/113.
assert 2*355*39894**2 < 113*100000**2
assert 2*333*39896**2 > 106*100000**2
glo,ghi=ratlo(39894,100000),rathi(39896,100000)
gL=[]; gU=[]
for j in range(XMAX*DX+1):
    x=j*S//DX
    ex=-(x*x//(2*S))
    gL.append(down(glo,elo(ex)))
    gU.append(up(ghi,ehi(ex)))
tail_g=gU[-1]
def density_upper_at_lower_arg(x):
    return up(ghi,ehi(-(x*x//(2*S))))
def monos(x,scale):
    el=[S]; eu=[S]
    fac_lo=elo(-x//scale); fac_hi=ehi(-x//scale)
    for _ in range(XMAX*scale):
        el.append(down(el[-1],fac_lo)); eu.append(up(eu[-1],fac_hi))
    return el,eu
def boxes():
    zl_lo,zl_hi=monos(0,DX)
    min_bound,max_arg=None,None
    panels=[]
    for i in range(ZMAX_N):
        z0=i*S//DZ; z1=(i+1)*S//DZ
        h0=S-down(z0,z0); h1=S-down(z1,z1)
        zr_lo,zr_hi=monos(z1,DX)
        aL=aU=mL=bU=0
        for j in range(DX*XMAX):
            x0=j*S//DX; x1=(j+1)*S//DX
            # V=(1-exp(-zx))/z+exp(-zx)*min(x/h,1/z).
            firstL=divlo(S-zr_hi[j],z1)
            firstU=divhi(S-zl_lo[j+1],z0) if z0 else x1
            secondL=down(zr_lo[j+1],min(divlo(x0,h0),divlo(S,z1)))
            secondU=divhi(x1,h1)
            if z0: secondU=min(secondU,divhi(S,z0))
            secondU=up(zl_hi[j],secondU)
            vL=max(0,firstL+secondL); vU=firstU+secondU
            # Sums before the common mesh-width scaling.
            aL+=down(vL,gL[j+1]); aU+=up(vU,gU[j])
            mL+=down(down(x0,vL),gL[j+1])
            bU+=up(up(vU,vU),gU[j])
        aL//=DX; aU=ceildiv(aU,DX)+3*tail_g
        mL//=DX; bU=ceildiv(bU,DX)+9*up(rathi(145,12),tail_g)
        # V<=3x past 12; integral x²g=(12+Mills(12))*g(12).
        kl=trans("log",divhi(bU,down(aL,aL)),True)
        assert kl>0
        U=divlo(down(mL,mL),up(aU,kl))
        q0=down(z0,z0)//2; q1=down(z1,z1)//2
        lamL=down(h1,elo(-q1)); lamU=up(h0,ehi(-q0))
        base=2*down(glo,elo(-q1))//3
        ent=(down(lamL,U)-2*up(lamU,aU))//2
        # Loss <= z²/(2h)g(h/z) + z²/(1+z²)g(1/z).
        zsq=up(z1,z1)
        loss1=up(divhi(zsq,2*h1),density_upper_at_lower_arg(divlo(h1,z1)))
        loss2=up(divhi(zsq,S+zsq),density_upper_at_lower_arg(divlo(S,z1)))
        bound=base+ent-loss1-loss2
        if min_bound is None or bound<min_bound: min_bound,max_arg=bound,i
        if i%(DZ//16)==0: panels.append(bound)
        else: panels[-1]=min(panels[-1],bound)
        zl_lo,zl_hi=zr_lo,zr_hi
    return min_bound,max_arg,panels
def large_z():
    z=ratlo(9,16); total=0
    for j in range(DX*XMAX):
        x=j*S//DX; a=down(z,x)
        # exp(-a)(1+a) and 2exp(-(1+a)/2) join continuously,
        # both decrease. The threshold is exact on this dyadic mesh.
        envelope=2*ehi(-(S+a)//2) if a<S else up(S+a,ehi(-a))
        total+=up(envelope,gU[j])
    # g-tail <=g(12)/12 and envelope<=2.
    return 2*ceildiv(total,DX)+ceildiv(4*tail_g,12)
def main():
    b,i,panels=boxes()
    env=large_z()
    print(json.dumps(dict(python=platform.python_version(),decimal_digits=getcontext().prec,
        scale=S,xmesh=DX,zmesh=DZ,xmax=XMAX,zmax="9/16",cells=DX*XMAX*ZMAX_N,
        minimum_lower_numerator=b,minimum_at_zindex=i,
        block_lower_numerators=panels,large_z_envelope_upper_numerator=env,
        success=b>0 and env<S),indent=2))
    assert b>0 and env<S
if __name__=="__main__": main()
