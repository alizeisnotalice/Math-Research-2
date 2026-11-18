"""Additional integer witnesses certified inside bodies contracted by 0.999999006."""
from fractions import Fraction as F
import numpy as np


def scaled_material(base,lower,upper,original):
    """Keep every old point and add certified points on interleaved angular samples."""
    from flint import arb,ctx
    ctx.prec=160
    scale=arb(999999006)/1000000000
    q=2**28
    rows=1025
    flat,counts=base
    assert list(map(int,counts))==[96,192,320,448,576]
    counts2=np.asarray([1024,384,640,896,1024],dtype=np.int32)
    assert max(counts2)<=1024 and sum(counts2)<=4096
    old=flat.reshape(rows,-1,2)
    dense=np.zeros((rows,sum(counts2),2),dtype=np.int32)
    cursor=offset=0
    pi=arb.pi()
    for b,n in enumerate((0,3,5,7,9)):
        total=int(counts[b])
        added=int(counts2[b])-total
        dense[:,cursor:cursor+total]=old[:,offset:offset+total]
        rad=arb(1)/2 if not n else 1/(2*(pi/(2*n)).cos())
        corners=[] if not n else [(rad*(2*pi*j/n).cos(),rad*(2*pi*j/n).sin()) for j in range(n)]
        assert not n or rad<1  # All unfilled rows use the scaled body's origin.
        extra=[]
        for k in range(added):
            theta=F(2*k+1,2*added) if not n else F(2*(k*total//added)+1,2*total)
            z=2*pi*theta.numerator/theta.denominator
            u,v=z.cos(),z.sin()
            if n:
                corner=(theta*n+F(n,2)+F(1,2)).__floor__()%n
                x,y=corners[corner]
                dot=u*x+v*y
                radius=dot+(1-rad*rad+dot*dot).sqrt()
            else:
                radius=rad
            shrink=1-arb(3)/10**6
            extra.append((u*radius*shrink,v*radius*shrink))
        if b<2:
            angles=range(1)
            hi_angle=0
        else:
            k=3*b-4
            top=int(original[k])
            low=max(0,(1024*int(lower[k]))//top-1)
            high=min(1024,-((-1024*int(upper[k]))//top)+1)
            angles=range(low,high+1)
            hi_angle=arb(top)/(2**44)
        for row in angles:
            theta=hi_angle*row/1024 if b>=2 else arb(0)
            c,s=theta.cos(),theta.sin()
            vertices=[(c*x-s*y,s*x+c*y)for x,y in corners]
            for j,(x,y) in enumerate(extra):
                xx,yy=c*x-s*y,s*x+c*y
                ix,iy=round(float(xx.mid())*q),round(float(yy.mid())*q)
                X,Y=arb(ix)/q,arb(iy)/q
                if n:
                    assert all((X-scale*vx)**2+(Y-scale*vy)**2<scale*scale for vx,vy in vertices)
                else:
                    assert X*X+Y*Y<scale*scale/4
                dense[row,cursor+total+j]=(ix,iy)
        # Outside the selected angular rows, each added point is the common origin.
        if b<2:
            dense[1:,cursor+total:cursor+total+added]=dense[0,cursor+total:cursor+total+added]
        cursor+=total+added
        offset+=total
    assert dense.shape==(1025,3968,2) and np.max(np.abs(dense))<q
    return dense.ravel(),counts2
