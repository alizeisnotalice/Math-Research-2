"""A sliding set of diameter contacts in two independently released oblique caps."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from flint import arb, ctx, __version__ as flint_version


def run(comp, source, max_height):
    sys.path.insert(0,str(comp.resolve()));sys.path.insert(0,str(source.resolve()))
    from certify_compulsory_C3 import EXTRA, UNIT
    import audit_optimized_eh_chord as g
    from certify_three_sigma_area_budget import pack, exact
    from certify_sigma_mandatory_arcs import degree
    ctx.prec=192;g.PI=arb.pi()
    half=arb(1)/2;face_slack=exact('.0058')
    def test(sl,sr,bl,br):
        g.SIGMA=degree(sl).union(degree(sr));b=exact(bl).union(exact(br))
        n,m,*_=g.setup();base=g.line_pair(n[0],m[5]);c=g.plus(base,(arb(0),b))
        v=g.plus(g.scale(half,n[1]),g.scale(exact('-.04'),g.rot(n[1])))
        pts={'C':c,'L':g.minus(c,m[5]),'V':v,'W':g.minus(v,n[1]),'A':(arb(0),arb('.501'))}
        pts.update({f'k{i}':tuple(exact(z) for z in p) for i,p in enumerate(EXTRA)})
        pts.update({'J'+k:(-p[0],p[1]) for k,p in list(pts.items()) if k!='A'})
        d2=[(1-g.dot(g.minus(p,q),g.minus(p,q))).lower() for (i,p),(j,q) in combinations(pts.items(),2) if frozenset((i,j)) not in UNIT]
        assert len(d2)==205 and min(d2)>0,('split',min(d2))
        def dirs(l,r):return [g.unit(g.PI+degree(l).union(degree(r))+2*j*g.PI/3) for j in range(3)]
        wind=sum(g.dot(g.rot(nn),pts[k]) for nn,k in zip(dirs(F(0),F(0)),('JC','JW','V')))
        end_slope=wind*degree(1).cos()-3*half*degree(1).sin()
        assert end_slope>0
        cap=min(max(g.dot(m[j],p).lower() for p in pts.values())-half for j in (0,1))
        assert cap>face_slack
        trig=-m[5][1]
        assert trig>exact('.4765') and trig<exact('.478')
        todo=[(F(1),F(30))];rooms=[]
        while todo:
            l,r=todo.pop()
            v=sum(max(g.dot(n,p).lower() for p in pts.values()) for n in dirs(l,r))-3*half
            if v>0:rooms.append(v.lower());continue
            assert r-l>F(1,100000),('split',l,r)
            mid=(l+r)/2;todo.extend(((l,mid),(mid,r)))
        return dict(sigma=list(map(str,(sl,sr))),height=list(map(str,(bl,br))),root_cells=len(rooms),
                    slacks={k:pack(v)['exact'][0] for k,v in
                            [('diameter',min(d2)),('root',min(rooms)),('endpoint',end_slope.lower()),('strong_cap',cap.lower()),('sine',trig.lower())]})
    def cover(sl,sr,bl,br):
        try:data.append(test(sl,sr,bl,br));return
        except AssertionError as e:
            if sr-sl<F('0.000025') and br-bl<F('.0000001'):raise
        if (sr-sl)/F('.04') > (br-bl)/max_height/12:
            mid=(sl+sr)/2;cover(sl,mid,bl,br);cover(mid,sr,bl,br)
        else:
            mid=(bl+br)/2;cover(sl,sr,bl,mid);cover(sl,sr,mid,br)
    data=[];cover(F('1.50'),F('1.54'),F(0),max_height)
    # The subdivision consists of exact binary partitions; check every vertical cross-section again.
    ticks=sorted({F(x) for r in data for x in r['sigma']})
    assert (ticks[0],ticks[-1])==(F('1.50'),F('1.54'))
    for l,r in zip(ticks,ticks[1:]):
        mids=(l+r)/2
        ints=sorted(tuple(map(F,z['height'])) for z in data if F(z['sigma'][0])<=mids<=F(z['sigma'][1]))
        assert ints[0][0]==0 and ints[-1][1]==max_height and all(a[1]==b[0] for a,b in zip(ints,ints[1:]))
    bounds={k:str(min(F(r['slacks'][k]) for r in data)) for k in data[0]['slacks']}
    return dict(sigma_degrees=['1.50','1.54'],height=[str(0),str(max_height)],minimum=bounds,
        angle_cells=sum(r['root_cells'] for r in data),parameter_rectangles=len(data),cells=data,
        positive_cap_indices=[0,1,2],checked_unit_pairs=[sorted(p) for p in sorted(UNIT,key=lambda z:sorted(z))],
        parent_sha256=sha256((source/'certify_compulsory_C3.py').read_bytes()).hexdigest(),
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),python_flint=flint_version,arb_bits=ctx.prec)
if __name__=='__main__':
    ap=ArgumentParser()
    for k in ('comp','prior','output'):ap.add_argument(k,type=Path)
    ap.add_argument('--height',type=F,default=F('.00011'))
    a=ap.parse_args();v=run(a.comp,a.prior,a.height);a.output.write_text(json.dumps(v,indent=1)+'\n')
    print(v['parameter_rectangles'],v['angle_cells'],{k:float(F(x)) for k,x in v['minimum'].items()})
