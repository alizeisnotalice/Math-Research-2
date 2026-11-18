import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys


def build(comp, seven, k=512, five_only=False, omit=(), fraction=None):
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
    p=list(map(lambda x:ar(F(x)),pose))
    pts=[]
    pi=arb.pi()
    for i in range(8*k):
        a=pi*i/(4*k);append_point(pts,a.cos()/2,a.sin()/2)
    if 3 not in omit: arc_body(pts,k,pi,3,arb(0),*p[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        if n not in omit: arc_body(pts,k,pi,n,*p[j:j+3])
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
