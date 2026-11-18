"""Generate R130's exact joint-cap kernel with piecewise angular support bounds."""
from fractions import Fraction as F
from pathlib import Path

import numpy as np

DIV=64


def piecewise_support(proof,box,Q,U,ceil):
    groups=[]
    for n,k in ((3,None),(5,2),(7,5),(9,8)):
        a,b=(F(0),F(0)) if k is None else box[k]
        step=(b-a)/DIV
        pad=ceil(Q*F(3,5)*step**2*F(U+2,U)/8)
        groups.extend(proof.one(n,a+i*step)+pad for i in range(DIV+1))
    return proof.directions,proof.disk,np.ascontiguousarray(groups,dtype=np.int64).ravel(),proof.pointers,proof.triples


def write_kernel(work,source):
    header=work/'fan-exposed.hpp'
    text=header.read_text()
    old='hi[index]>lo[index] && (hi[index]+lo[index])%2==0'
    assert text.count(old)==1
    (work/'fan-mesh.hpp').write_text(text.replace(old,'hi[index]-lo[index]>=2'))
    code=source.read_text().replace('"fan-exposed.hpp"','"fan-mesh.hpp"')
    before='''        auto support_upper=[&](int j){'''
    grid='''        constexpr int subdivisions=64;
        int ix0[4]{},ix1[4]{};
        L frac0[4]{},frac1[4]{},den[4]{};
        for(int i=1;i<4;++i){
            int k=3*i-1;
            den[i]=root_hi[k]-root_lo[k];
            if(den[i]<=0 || lo[k]<root_lo[k] || hi[k]>root_hi[k])std::abort();
            L p0=(lo[k]-root_lo[k])*subdivisions,p1=(hi[k]-root_lo[k])*subdivisions;
            ix0[i]=std::min<L>(subdivisions-1,p0/den[i]);
            ix1[i]=std::min<L>(subdivisions-1,p1/den[i]);
            frac0[i]=p0-ix0[i]*den[i];frac1[i]=p1-ix1[i]*den[i];
        }
''' + before
    assert code.count(before)==1
    code=code.replace(before,grid)
    start=code.index('                L h=endpoints[(2*i)*n+j],p=endpoints[(2*i+1)*n+j];')
    end=code.index('                int k=3*i;',start)
    code=code[:start]+'''                const L* values=endpoints+i*(subdivisions+1)*n+j;
                L h=values[0];
                if(i){
                    auto estimate=[&](int s,L fraction){
                        return up_signed(W(values[s*n])*(den[i]-fraction)+W(values[(s+1)*n])*fraction,den[i]);
                    };
                    h=std::max(estimate(ix0[i],frac0[i]),estimate(ix1[i],frac1[i]));
                    for(int s=ix0[i]+1;s<=ix1[i];++s)h=std::max(h,values[s*n]);
                }
'''+code[end:]
    before='        if(b<=a || (a+b)%2)return stop(lo,hi,3);'
    assert code.count(before)==1
    code=code.replace(before,'        if(b-a<2)return stop(lo,hi,3);')
    result=work/'weighted-mesh.cpp'
    result.write_text(code)
    return result


def dense_material(base,lower,upper,original):
    """Keep every old point and add certified points on interleaved angular samples."""
    from flint import arb,ctx
    ctx.prec=160
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
            shrink=1-arb(1)/10**7
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
                    assert all((X-vx)**2+(Y-vy)**2<1 for vx,vy in vertices)
                else:
                    assert X*X+Y*Y<arb(1)/4
                dense[row,cursor+total+j]=(ix,iy)
        # Outside the selected angular rows, each added point is the common origin.
        if b<2:
            dense[1:,cursor+total:cursor+total+added]=dense[0,cursor+total:cursor+total+added]
        cursor+=total+added
        offset+=total
    assert dense.shape==(1025,3968,2) and np.max(np.abs(dense))<q
    return dense.ravel(),counts2
