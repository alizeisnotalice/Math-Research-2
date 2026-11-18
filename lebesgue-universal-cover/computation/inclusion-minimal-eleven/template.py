import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys


def build(comp, seven, k=512, five_only=False, omit=(), fraction=None, tube_n=None, extra=0, box=None):
    assert sha256((Path(comp)/'certify_nine.py').read_bytes()).hexdigest()=='610e94cf22e901453bc34f77afb8f0994d9e471432cbe20498fca3b1530c9541'
    sys.path.insert(0,str(Path(comp).resolve()))
    from certify_nine import S, append_point, ar, arb, arc_body, ctx, hull
    ctx.prec=160
    def read(name):
        pins={
            'free-p30-reoptimized-pose.json':'21e6fceeb7f2033b69a6c2a7ae633abb26b6afd3621e32f6e6aa57f9b8770593',
            'seven-point-families-pose.json':'114b8e680bb0adaeb6453b1a87abe856616873e32dd6c01a946b3291e917ee8b',
            'finite-point-leaf-input.json':'b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc',
            'multi-host-points.json':'9052ca388aabacc87116ba99a05d9ef98b49dfeed564e8cfd18efb647f93f282',
        }
        raw=(Path(comp)/name).read_bytes()
        assert sha256(raw).hexdigest()==pins[name]
        return json.loads(raw)
    pose=read('seven-point-families-pose.json')['pose'] if seven else read('free-p30-reoptimized-pose.json')
    if fraction is not None:
        assert five_only and not omit
        low=read('free-p30-reoptimized-pose.json')
        high=read('seven-point-families-pose.json')['pose']
        t=F(fraction)
        pose=[(1-t)*F(x)+t*F(y) for x,y in zip(low[:11],high[:11])]
    if box is not None:
        assert len(box)==11 and all(type(l) is int and type(h) is int and l<=h for l,h in box)
        assert five_only and not omit and fraction is None and tube_n is None
        pose=[F(l+h,2**45) for l,h in box]
    p=list(map(lambda x:ar(F(x)),pose))
    # A moving B_n contains the midpoint copy contracted about its midpoint center
    # by displacement divided by the inradius. The disk itself remains fixed.
    steps={}
    if box is not None:
        half=[F(h-l,2**45) for l,h in box]
        for n,j in ((3,None),(5,2),(7,5),(9,8)):
            pos=0 if j is None else j+1
            motion=ar(half[pos]**2+half[pos+1]**2).sqrt()
            if j is not None:motion+=(arb.pi()/(2*n)).tan()/2*ar(half[j])
            steps[n]=1-motion/(1-1/(2*(arb.pi()/(2*n)).cos()))
            assert 0<steps[n]<=1
    if tube_n:
        assert fraction is not None and five_only and not omit and F(extra)>=0
        for n,j in ((3,None),(5,2),(7,5),(9,8)):
            pos=0 if j is None else j+1
            dx,dy=(F(high[pos+i])-F(low[pos+i]) for i in range(2))
            speed=ar(dx*dx+dy*dy).sqrt()
            if j is not None:
                speed+=(arb.pi()/(2*n)).tan()/2*ar(abs(F(high[j])-F(low[j])))
            motion=speed/(2*tube_n)+ar(F(extra))
            radius=1-1/(2*(arb.pi()/(2*n)).cos())
            steps[n]=1-motion/radius
            assert 0<steps[n]<1
    pts=[]
    pi=arb.pi()
    for i in range(8*k):
        a=pi*i/(4*k);append_point(pts,a.cos()/2,a.sin()/2)
    def add_body(n,rho,tx,ty):
        if n not in steps:
            arc_body(pts,k,pi,n,rho,tx,ty)
            return
        scale=steps[n]
        radius=scale/(2*(pi/(2*n)).cos())
        for ell in range(n):
            phi=rho+2*ell*pi/n
            cx,cy=tx+radius*phi.cos(),ty+radius*phi.sin()
            for j in range(k+1):
                ang=phi+pi*(1-ar(F(1,2*n))+ar(F(j,n*k)))
                append_point(pts,cx+scale*ang.cos(),cy+scale*ang.sin())
    if 3 not in omit:add_body(3,arb(0),*p[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        if n not in omit:add_body(n,*p[j:j+3])
    for name,den,key in [('finite-point-leaf-input.json','denominator',11),('multi-host-points.json','point_denominator',14)]:
        if five_only or (key==14 and not seven): continue
        cloud=read(name)
        xy=[(F(x,cloud[den]),F(y,cloud[den])) for x,y in cloud['points']]
        if key==14:
            a,b=[sum(row[j] for row in xy)/len(xy) for j in range(2)]
            xy=[(x-a,(y-b)*F(pose[17])) for x,y in xy]
        sn,cs=p[key].sin_cos()
        for x,y in xy:append_point(pts,cs*ar(x)-sn*ar(y)+p[key+1],sn*ar(x)+cs*ar(y)+p[key+2])
    return hull(pts),S+4
