"""A vertical-symmetric diameter-one anchor family forces the Gibbs corner C3 in every universal subhost of P(sigma)."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json,sys
from flint import arb,ctx
EXTRA=(
 ('.026049970901','.515888228177'),('.336263316958','-.347577565759'),
 ('.181470788663','.464793989670'),('.447623842314','-.212397700613'),
 ('.293446043323','-.390014410874'),('.067495484760','.504591212874'))
UNIT={frozenset(x) for x in [('C','JC'),('C','L'),('JC','JL'),('V','W'),('JV','JW')]}
def run(comp):
    sys.path.insert(0,str(comp))
    import audit_optimized_eh_chord as g
    from certify_three_sigma_area_budget import pack,exact
    from certify_sigma_mandatory_arcs import degree
    ctx.prec=192
    g.PI=arb.pi()
    lam=exact(F(1,10**7));half=arb('.5')
    def test(a,b):
        g.SIGMA=degree(a).union(degree(b))
        n,m,*_=g.setup();C=g.line_pair(n[0],m[5]);L=g.minus(C,m[5])
        V=g.plus(g.scale(half,n[1]),g.scale(exact('-.04'),g.rot(n[1])));W=g.minus(V,n[1])
        pts={'C':C,'L':L,'V':V,'W':W,'A':(arb(0),arb('.501'))}
        pts.update({f'k{i}':tuple(exact(v) for v in xy) for i,xy in enumerate(EXTRA)})
        pts.update({'J'+k:(-v[0],v[1]) for k,v in list(pts.items()) if k!='A'})
        assert len(pts)==21
        diam=[(1-g.dot(g.minus(v,z),g.minus(v,z))).lower()
              for (k,v),(l,z) in combinations(pts.items(),2) if frozenset((k,l)) not in UNIT]
        assert len(diam)==205 and min(diam)>0,('diameter',a,b,min(diam))
        # L=C-m5 gives m2.L exactly 1/2. A positive harmonic adds cos(3sigma)/100 there.
        cap=[(1-lam)*(max(g.dot(m[i],p).lower() for p in pts.values())-half)-lam/100 for i in (0,1)]
        assert min(cap)>0 and (3*g.SIGMA).cos()>0
        def dirs(x,y):return [g.unit(g.PI+degree(x).union(degree(y))+2*i*g.PI/3) for i in range(3)]
        at_zero=dirs(F(0),F(0))
        B=sum(g.dot(g.rot(v),pts[k]) for v,k in zip(at_zero,('JC','JW','V'),strict=True))
        one=degree(F(1))
        derivative=(1-lam)*(B*one.cos()-3*half*one.sin())-9*lam/100
        assert derivative>0
        passed=[]
        def scan(x,y):
            dirs0=dirs(x,y)
            lower=sum(max(g.dot(u,p).lower() for p in pts.values()) for u in dirs0)-3*half
            room=(1-lam)*lower-3*lam/100
            if room>0:
                passed.append(room.lower());return
            assert y-x>F(1,1000000),(a,b,x,y,room)
            mid=(x+y)/2;scan(x,mid);scan(mid,y)
        scan(F(1),F(30))
        return dict(sigma=[str(a),str(b)],nonunit_diameter_slack=pack(min(diam)),
                    m0_m1_slack=pack(min(cap).lower()),m2_per_weight=pack(((3*g.SIGMA).cos()/100).lower()),
                    root_derivative=pack(derivative.lower()),root_boxes=len(passed),root_room=pack(min(passed)))
    data=[]
    def visit(a,b):
        try:r=test(a,b)
        except AssertionError:
            assert b-a>F(1,100000), (a,b)
            mid=(a+b)/2;visit(a,mid);visit(mid,b)
        else:data.append(r)
    visit(F('1.50'),F('1.54'))
    assert data and data[0]['sigma'][0]==str(F('1.50')) and data[-1]['sigma'][1]==str(F('1.54'))
    for k,l in zip(data,data[1:]):assert k['sigma'][1]==l['sigma'][0]
    wall=exact('.000000325193131685913632299487350053');c3=exact('.32512346')
    assert wall>exact('.00000100021')*c3
    return dict(complete=True,interval=['1.50','1.54'],input_pairs=EXTRA,beta='-.04',axis_height='.501',unit_pairs=sorted(sorted(p) for p in UNIT),
                total_angle_boxes=sum(r['root_boxes'] for r in data),sigma_bins=data,precision=ctx.prec,
                SHA256=sha256(Path(__file__).read_bytes()).hexdigest(),python=sys.version)

if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('comp',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result=run(a.comp.resolve());a.output.write_text(json.dumps(result,indent=1)+'\n')
    print(len(result['sigma_bins']),result['total_angle_boxes'])
