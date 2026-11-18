import argparse,sys,json
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from hashlib import sha256
sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build
from flint import arb,ctx

def q(z):
    z=F(z);return arb(z.numerator)/z.denominator

def run(comp,input,k,cloud=None):
    sys.path.insert(0,str(comp.resolve()));from certify_nine import append_point,hull
    ctx.prec=160
    row=json.loads(input.read_bytes());box=[[v,v] for v in row['center']]
    assert len(box)==11 and all(type(v[0]) is int for v in box)
    hh,D=build(comp,False,k,True,box=box);S=D-4
    point_source=cloud or (comp/'finite-point-leaf-input.json')
    raw=point_source.read_bytes()
    assert sha256(raw).hexdigest() in ('b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc','45cb554f79d2c30f23cdc526f7cc1a333b6f98c3dbaa60234862d449bb66de05')
    sources=[(json.loads(raw),row['point_pose'])]
    if 'second_point_pose' in row:
        assert cloud is None
        second=(Path(__file__).resolve().parent/'new-points.json').read_bytes()
        assert sha256(second).hexdigest()=='45cb554f79d2c30f23cdc526f7cc1a333b6f98c3dbaa60234862d449bb66de05'
        sources.append((json.loads(second),row['second_point_pose']))
    pts=list(hh)
    for pointset,pose in sources:
        den=pointset['denominator']
        x,y,a=map(q,pose);s,c=a.sin_cos()
        for X,Y in pointset['points']:
            X=q(F(X,den));Y=q(F(Y,den));append_point(pts,x+c*X-s*Y,y+s*X+c*Y)
    pg=hull(pts);lo=F(sum(p[0]*r[1]-p[1]*r[0] for p,r in zip(pg,pg[1:]+pg[:1])),2*D**2)
    eps=F(2,S)+F(22,7)**2/(72*k*k)
    per=sum(1+isqrt((p[0]-r[0])**2+(p[1]-r[1])**2) for p,r in zip(pg,pg[1:]+pg[:1]))
    hi=lo*F(D,S)**2+eps*F(per,S)+F(22,7)*eps**2
    return dict(lower=str(lo),upper=str(hi),decimal=[f'{float(lo):.15f}',f'{float(hi):.15f}'],k=k,pose_sha256=sha256(input.read_bytes()).hexdigest())
if __name__=='__main__':
 p=argparse.ArgumentParser()
 for t in ('comp','pose','output'):p.add_argument(t,type=Path)
 p.add_argument('--k',type=int,default=4096);p.add_argument('--cloud',type=Path);args=p.parse_args();d=run(args.comp,args.pose,args.k,args.cloud)
 args.output.write_text(json.dumps(d,indent=2)+'\n');print(d['decimal'])
