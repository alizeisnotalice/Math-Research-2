"""Certify five-point containment in the moving five-input host and a directed area upper bound."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build
from flint import arb,ctx


def q(value):
    z=F(value);return arb(z.numerator)/z.denominator


def verify(comp,points,center,pose,hand,k):
    assert hand in [-1,1] and len(center)==11 and all(type(x) is int for x in center)
    all_bounds=[[c,c] for c in center]
    h,D=build(comp,False,k,True,box=all_bounds)
    S=D-4
    ctx.prec=160
    x,y,angle=map(q,pose)
    sn,cs=angle.sin_cos()
    five=json.loads(points.read_bytes())
    assert five['indices']==[0,2,5,7,9] and len(five['points'])==5
    p=[(q(a),q(b)*hand) for a,b in five['points']]
    rotated=[((x+cs*a-sn*b)*D,(y+sn*a+cs*b)*D) for a,b in p]
    verts=h[1:]+h[:1]
    dets=[];twice=per=0
    for (a,b),(c,d) in zip(h,verts):
        dx=c-a;dy=d-b;twice+=a*d-b*c;per+=1+isqrt(dx*dx+dy*dy)
        dets.append(min((dx*(Y-b)-dy*(X-a)).lower() for X,Y in rotated))
    margin=min(dets)/(D*D)
    eps=F(2,S)+F(22,7)**2/F(72*k*k)
    lower=F(twice,2*D*D)
    upper=F(twice,2*S*S)+eps*F(per,S)+F(22,7)*eps*eps
    return margin,lower,upper


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('comp','points','input','out'):p.add_argument(key,type=Path)
    p.add_argument('--one',action='store_true');p.add_argument('--k',type=int,default=512)
    a=p.parse_args();assert __debug__ and a.k>=64
    data=json.loads(a.input.read_bytes())
    out=[]
    for row in data['rows']:
        box=row['bounds'];center=row['center']
        assert all(l<=c<=u for (l,u),c in zip(box,center,strict=True))
        assert all(l<=x<=y<=u for (l,u),(x,y) in zip(data['original_bounds'],box,strict=True))
        m,lo,hi=verify(a.comp,a.points,center,row['pose'],row['hand'],a.k)
        out.append(dict(label=row['label'],contained=bool(m>0),directed_edge_determinant_lower=str(m),
                        host_lower=str(lo),host_upper=str(hi),upper_decimal=str(F((hi*10**12).__ceil__(),10**12)),
                        below_083750=bool(m>0 and hi<F('.83750'))))
        if a.one:break
    summary=dict(arc_samples=a.k,rows=out,certified_feasible_below_083750=sum(r['below_083750'] for r in out),
                 original_leaf=data['leaf'],source_sha256=sha256(a.input.read_bytes()).hexdigest(),
                 certificate_role='explicit compatible placement inside each indicated original subbox',
                 points_sha256=sha256(a.points.read_bytes()).hexdigest())
    a.out.write_text(json.dumps(summary,indent=2)+'\n')
    print(summary['certified_feasible_below_083750'],'/',len(out))
    for r in sorted((r for r in out if r['contained']),key=lambda r:F(r['host_upper']))[:2]:print(r['label'],float(F(r['host_upper'])),r['upper_decimal'],flush=True)
if __name__=='__main__':main()
