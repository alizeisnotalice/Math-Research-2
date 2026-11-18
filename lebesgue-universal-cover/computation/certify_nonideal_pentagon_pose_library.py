"""Enclose a whole parameter box of genuine pentagons and exclude all but two poses in P."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from flint import arb, ctx
from deformed_reuleaux_support import body, environment


def verify(ns, w, e, rw, re):
    A, box, mul, dot, perp = (ns[k] for k in ('A','box','mul','dot','perp'))
    vertices, support = body(ns,box(w-rw,w+rw),box(e-re,e+re))
    edges = [ns['sub'](vertices[i],vertices[j]) for i,j in ns['orders']]
    cross = ns['cross']
    normals = list(map(ns['av'],ns['d']['n']))
    cap_normals = list(map(ns['av'],ns['d']['m']))
    pins = []
    for k,i in enumerate((1,3,0)):
        u = mul(-1,normals[k])
        slots = [j for j,v in enumerate(edges) if cross(v,u)>0 and cross(u,edges[(j+1)%10])>0]
        assert len(slots)==1 and tuple(ns['sector'][slots[0]])==(i,0)
        assert dot(vertices[i],u).contains(A('.5'))
        pins.append(slots[0])
    slope = sum((-1)**k*dot(vertices[i],perp(normals[k])) for k,i in enumerate((1,3,0)))
    sign = 1 if slope>0 else -1
    endpoint = F('.01')
    endpoint_margin = (sign*slope-A('1.5')*(A(endpoint/2)*ns['angle_unit']).tan()).lower()
    assert endpoint_margin>0, 'endpoint slope'
    minimum, count, max_depth = A(1), 0, 0
    def visit(l,r,depth):
        nonlocal minimum,count,max_depth
        mid = (l+r)/2
        fm = ns['F_at'](vertices,edges,mid)
        assert not sign*fm<0, 'interior sign reversal'
        derivative = sum((-1)**k*ns['supports'](vertices,edges,ns['trig'](box(l+60*k,r+60*k)),normals[0],True) for k in range(3))
        value = (sign*(fm+box(l-mid,r-mid)*ns['angle_unit']*derivative)).lower()
        if value>0:
            minimum=min(minimum,value); count+=1; max_depth=max(max_depth,depth)
            return
        assert depth<14, 'inconclusive interior'
        visit(l,mid,depth+1)
        visit(mid,r,depth+1)
    for k in range(240):
        visit(endpoint+(60-2*endpoint)*k/240,endpoint+(60-2*endpoint)*(k+1)/240,0)
    losses=[]
    for turn in range(0,360,60):
        for s in (1,-1):
            if (turn,s) in ((0,1),(240,-1)):
                continue
            pair=[]
            for k in (3,5):
                u=cap_normals[k]
                height=ns['supports'](vertices,edges,ns['trig'](A(-turn)),(u[0],s*u[1]))
                pair.append(((height-A('.5')).lower(), k))
            gap,edge=max(pair)
            assert gap>0, 'cap'
            losses.append([turn,s,edge,str(gap)])
    return {'w':str(w),'deformation':str(e),'w_radius':str(rw),'e_radius':str(re),'pinned_sectors':pins,'F_sign':sign,'endpoint_margin':str(endpoint_margin),'interior_margin':str(minimum),'interior_cells':count,'max_depth':max_depth,'cap_eliminations':losses}


def run(archive,inp,out):
    ctx.prec=140
    cfg=json.loads(inp.read_bytes())
    ns=environment(str(archive.resolve()),cfg['sigma'])
    rows=[];failed=[]
    for ws,es in cfg['points']:
        w,e=F(ws),F(es)
        for rw,re in [(F(cfg['w_radius']),F(cfg['e_radius'])),(F(cfg['w_radius'])/10,F(cfg['e_radius'])/10),(F(0),F(0))]:
            try:
                row=verify(ns,w,e,rw,re)
                rows.append(row)
                break
            except (ArithmeticError,AssertionError) as exc:
                reason=str(exc)
        else:
            failed.append([ws,es,reason])
    sources=[Path(__file__),Path(__file__).with_name('deformed_reuleaux_support.py'),Path(__file__).with_name('certify_sigma_extended_graph.py')]
    obj={'parameters':cfg,'proofs':rows,'unresolved':failed,'scope':'For each closed box only id and Q_A can place the body into P.','source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in sources}}
    out.write_text(json.dumps(obj,indent=2)+'\n')
    print(len(rows),len(failed),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('parameters',type=Path);p.add_argument('out',type=Path);args=p.parse_args();run(args.archive,args.parameters,args.out)
