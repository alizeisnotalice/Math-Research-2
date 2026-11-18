"""Sample exact Arb support points of the 117 posed pentagons and output their integer convex hull."""
import argparse, hashlib, json, math, platform, sys, time
from collections import Counter
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
import flint

GRID=10**11

def hull(points):
    pts=sorted(set(points))
    def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def half(seq):
        out=[]
        for p in seq:
            while len(out)>=2 and cross(out[-2],out[-1],p)<=0:out.pop()
            out.append(p)
        return out
    if len(pts)<2:return pts
    return half(pts)[:-1]+half(reversed(pts))[:-1]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(computation,poses_file,out,N):
    assert N>=4
    computation=computation.resolve();sys.path.insert(0,str(computation))
    from deformed_reuleaux_support import body,environment
    ctx.prec=160
    records=json.loads(poses_file.read_bytes())['poses']
    expected={'nonideal-pentagon-library.json','nonmonotone-pentagon-inputs.json'}
    records=[x for x in records if x['source'] in expected]
    assert len(records)==117 and len({x['index'] for x in records})==117
    inputs={name:json.loads((computation/name).read_bytes()) for name in expected}
    for r in records:
        idx=r['body_index'];item=inputs[r['source']]
        if r['source']=='nonideal-pentagon-library.json':v=item['proofs'][idx];a,b=v['w'],v['deformation']
        else:a,b=item[idx]
        assert Fraction(a)==Fraction(r['w']) and Fraction(b)==Fraction(r['e'])
        assert r['reflection'] in (-1,1) and isinstance(r['turn'],int)
    ns=environment(str(computation),'1.5497');ctx.prec=160
    add,mul,cross=(ns[k] for k in ('add','mul','cross'))
    dirs=[]
    for k in range(N):
        sn,cs=(2*arb.pi()*k/N).sin_cos();dirs.append((cs,sn))
    trigonometry={}
    for key in {(x['turn'],x['reflection']) for x in records}:
        turn,sgn=key;sn,cs=(arb.pi()*turn/180).sin_cos()
        trigonometry[key]=(cs,sn,[(cs*x+sn*sgn*y,-sn*x+cs*sgn*y) for x,y in dirs])
    merged=[];summary=[];ambiguous=0;start=time.monotonic()
    for row in records:
        vertices,_=body(ns,Fraction(row['w']),Fraction(row['e']))
        edges=[ns['sub'](vertices[i],vertices[j]) for i,j in ns['orders']]
        ef=[(float(x.mid()),float(y.mid())) for x,y in edges]
        cs,sn,localdirs=trigonometry[(row['turn'],row['reflection'])];sgn=row['reflection']
        points=[];fallback=0
        def transform(v):return (cs*v[0]-sn*v[1],sgn*(sn*v[0]+cs*v[1]))
        for dx,dy in localdirs:
            x,y=float(dx.mid()),float(dy.mid())
            approx=[a*y-b*x for a,b in ef]
            idx=next((i for i in range(10) if approx[i]>0 and approx[(i+1)%10]<0),None)
            direction=(dx,dy)
            if idx is not None and cross(edges[idx],direction)>0 and cross(direction,edges[(idx+1)%10])>0:
                candidates=[idx]
            else:
                fallback+=1
                cc=[cross(e,direction) for e in edges]
                candidates=[i for i in range(10) if not cc[i]<0 and not cc[(i+1)%10]>0]
            assert candidates
            possible=[]
            for i in candidates:
                k,c=ns['sector'][i];possible.append(transform(add(vertices[k],mul(c,direction))))
            q=tuple(round(float(z.mid())*GRID) for z in possible[0])
            # All sector branches not rigorously excluded contain the true support branch.
            # Each coordinate of every possible point is within 1/GRID of q.
            for pt in possible:
                for v,n in zip(pt,q):assert abs(v*GRID-n)<1
            points.append(q)
        hp=hull(points);assert len(hp)>=3
        merged=hull(merged+hp);ambiguous+=fallback
        summary.append({'index':row['index'],'source':row['source'],'body_index':row['body_index'],'w':row['w'],'e':row['e'],'turn':row['turn'],'reflection':sgn,'samples':N,'distinct_samples':len(set(points)),'body_hull_vertices':len(hp),'ambiguous_sector_directions':fallback})
        if len(summary)%10==0:print('sampled',len(summary),'/',len(records),'merged_vertices',len(merged),'seconds',round(time.monotonic()-start,1),flush=True)
    packed=''.join(f'{x},{y}\n' for x,y in merged).encode()
    deps=['deformed_reuleaux_support.py','certify_sigma_extended_graph.py','certify_ah_two_line_cap.py','nonideal-pentagon-library.json','nonmonotone-pentagon-inputs.json']
    result={'schema':'posed-pentagons-integer-hull-v1','python':platform.python_version(),'python_flint':flint.__version__,'precision_bits':ctx.prec,'sigma_degrees':'1.5497','directions':N,'grid':GRID,'body_count':len(records),'sample_count':len(records)*N,'coordinate_rounding_error_strict':'1/100000000000','euclidean_rounding_error_upper':'1/50000000000','ambiguous_sector_directions':ambiguous,'transform_convention':'(x,y) -> (cos(k)x-sin(k)y, reflection*(sin(k)x+cos(k)y))','vertices_ccw_integer':merged,'vertices_sha256':hashlib.sha256(packed).hexdigest(),'per_body':summary,'poses_sha256':digest(poses_file),'source_sha256':{name:digest(computation/name) for name in deps},'script_sha256':digest(Path(__file__))}
    out.parent.mkdir(parents=True,exist_ok=True);tmp=out.with_suffix(out.suffix+'.partial');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)
    print('complete',len(records),'samples',len(records)*N,'hull',len(merged),'seconds',round(time.monotonic()-start,1),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('computation',type=Path);p.add_argument('poses',type=Path);p.add_argument('output',type=Path);p.add_argument('N',type=int,nargs='?',default=16384);a=p.parse_args();run(a.computation,a.poses,a.output,a.N)
