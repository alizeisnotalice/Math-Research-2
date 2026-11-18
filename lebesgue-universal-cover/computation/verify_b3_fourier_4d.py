#!/usr/bin/env python3
"""Independent exact verification of the reported 4D B3/Fourier fan leaves."""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys, time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

TAU=F(419,500)
def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def det(a,b): return a[0]*b[1]-a[1]*b[0]
def hull(points):
    tags={}
    for x,y,z in points: tags.setdefault((x,y),z)
    p=[(x,y,tags[(x,y)]) for x,y in sorted(tags)]
    def f(s):
        q=[]
        for x in s:
            while len(q)>1 and cross(q[-2],q[-1],x)<=0:q.pop()
            q.append(x)
        return q[:-1]
    return f(p)+f(reversed(p))
def b3q():
    q=[(F(11,20),F(0)),(F(-11,40),F(47,100)),(F(-11,40),F(-47,100))]
    lo,hi=F(1732050807568877,10**15),F(1732050807568878,10**15); assert lo*lo<3<hi*hi
    a,b=1/hi,1/lo
    c=[((a,b),(F(0),F(0))),((-b/2,-a/2),(F(1,2),F(1,2))),((-b/2,-a/2),(F(-1,2),F(-1,2)))]
    for x,y in q:
        for (xl,xh),(yl,yh) in c:
            assert max(abs(x-xl),abs(x-xh))**2+max(abs(y-yl),abs(y-yh))**2<1
    return q
def fq():
    r=F(283,600); assert F(1,2)-F(1,80)-F(19,1200)==r
    q=[(F(47,100),0),(F(33,100),F(33,100)),(0,F(47,100)),(-F(33,100),F(33,100)),(-F(47,100),0),(-F(33,100),-F(33,100)),(0,-F(47,100)),(F(33,100),-F(33,100))]
    assert all(x*x+y*y<r*r for x,y in q); return q
def main():
    ap=argparse.ArgumentParser();ap.add_argument("certificate",type=Path);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args(); started=time.monotonic()
    raw=a.certificate.read_bytes();r=json.loads(raw)
    domain=Path(__file__).with_name("outer_domain_certificate.json");dbytes=domain.read_bytes(); assert hashlib.sha256(dbytes).hexdigest()==r["source_global_box"]["sha256"]
    inp=json.loads(dbytes);den=int(r["coordinate_denominator"]);assert den==int(inp["disk_polygon_denominator"])
    disk=[(int(x),int(y),"D") for x,y in inp["disk_polygon_points"]];assert all(x*x+y*y<=den*den//4 for x,y,_ in disk)
    bs=[(int(x*den),int(y*den),f"B{i}") for i,(x,y) in enumerate(b3q())]; fs=[(int(x*den),int(y*den),f"F{i}") for i,(x,y) in enumerate(fq())]
    bb={z:(x,y) for x,y,z in bs};fb={z:(x,y) for x,y,z in fs};root=((-21*den//100,21*den//100),(-21*den//100,21*den//100),(-6*den//25,6*den//25),(-6*den//25,6*den//25));threshold=(2*TAU*den*den).numerator
    def boxof(path):
        box=list(root)
        for depth,ch in enumerate(path):
            u,v=box[depth%4];m=(u+v)//2;box[depth%4]=(u,m) if ch=="0" else (m,v)
        return tuple(box)
    def point(v,t):
        x,y,z=v
        if z=="D":return x,y
        if z.startswith("B"):return bb[z][0]+t[0],bb[z][1]+t[1]
        return fb[z][0]+t[2],fb[z][1]+t[3]
    leaves={z["path"]:z for z in r["leaves"]};assert len(leaves)==len(r["leaves"])
    maxd=int(r["max_depth"])
    def complete(p):
        if p in leaves:return True
        return len(p)<maxd and complete(p+"0") and complete(p+"1")
    assert complete("")
    # No terminal leaf may have a terminal descendant.
    for p in leaves:
        assert not any(q!=p and q.startswith(p) for q in leaves)
    # Recreate B3-only grid to independently verify the reported comparison labels.
    bp=Path(__file__).with_name("027-b3-pruning.json");braw=bp.read_bytes();assert hashlib.sha256(braw).hexdigest()==r["comparison_b3_certificate"]["sha256"]
    br=json.loads(braw);assert int(br["coordinate_denominator"])==den; assert tuple(map(int,br["translation_domain_scaled"]))==root[0]
    assert F(br["tau"])==TAU
    assert br["source_global_box"]["sha256"]==hashlib.sha256(dbytes).hexdigest()
    n=2**int(br["max_level"]);lo,hi=root[0];fine=(hi-lo)//n;g=bytearray(n*n)
    for z in br["leaves"]:
        if z["status"]=="prune":
            x0,x1=map(int,z["x_scaled"]);y0,y1=map(int,z["y_scaled"])
            for i in range((x0-lo)//fine,(x1-lo)//fine):
                for j in range((y0-lo)//fine,(y1-lo)//fine):g[i*n+j]=1
    def relation(box):
        (x0,x1),(y0,y1),_,_=box;hits=sum(g[i*n+j] for i in range((x0-lo)//fine,(x1-lo)//fine) for j in range((y0-lo)//fine,(y1-lo)//fine));total=((x1-x0)//fine)*((y1-y0)//fine)
        return "fully_b3_projection_pruned" if hits==total else ("disjoint_from_b3_projection" if hits==0 else "mixed_with_b3_projection")
    checked=0
    for path,z in leaves.items():
        box=boxof(path); assert [[str(u),str(v)] for u,v in box]==z["box_scaled"]
        if z["status"]!="prune":continue
        assert z["reason"]=="strict_inner_fan_area_exceeds_tau" and z["b3_projection_relation"]==relation(box)
        t0=tuple((u+v)//2 for u,v in box);h=hull(disk+[(x+t0[0],y+t0[1],tag) for x,y,tag in bs]+[(x+t0[2],y+t0[3],tag) for x,y,tag in fs])
        assert all(point(v,t0)==(v[0],v[1]) for v in h);assert len(h)==z["fan_vertices_at_center"]
        vals=[]
        for t in itertools.product(*box):
            pts=[point(v,t) for v in h];assert all(det(x,y)>0 for x,y in zip(pts,pts[1:]+pts[:1]));vals.append(sum(det(x,y) for x,y in zip(pts,pts[1:]+pts[:1])))
        assert min(vals)==int(z["minimum_twice_fan_area_scaled"])>threshold;checked+=1
    assert all(z["status"] in ("prune","retain") for z in leaves.values())
    pruned=[z for z in leaves.values() if z["status"]=="prune"]
    measure=lambda rows:sum((F(1,2**len(z["path"])) for z in rows),F(0))
    pv=measure(pruned); rv=measure([z for z in leaves.values() if z["status"]=="retain"])
    assert pv+rv==1
    assert r["terminal_leaf_counts"]=={"all":len(leaves),"pruned":len(pruned),"retained":len(leaves)-len(pruned)}
    vf=r["four_dimensional_volume_fraction"]
    assert F(vf["pruned"])==pv and F(vf["retained"])==rv
    cnt=Counter(z["b3_projection_relation"] for z in pruned)
    assert vf["pruned_by_relation_to_B3_projection_counts"]==dict(cnt)
    relvol={k:str(measure([z for z in pruned if z["b3_projection_relation"]==k])) for k in cnt}
    assert vf["pruned_by_relation_to_B3_projection_volume_fraction"]==relvol
    out={"status":"verified","certificate":str(a.certificate),"certificate_sha256":hashlib.sha256(raw).hexdigest(),"checked_pruned_leaves":checked,"recomputed_pruned_fraction":str(pv),"recomputed_retained_fraction":str(rv),"recomputed_comparison_fractions":relvol,"terminal_leaf_count":len(leaves),"partition":"exact complete binary prefix tree from the 4D root box","scope":"checks reported 4D exclusions only; all retained leaves and other F8 coordinates remain unresolved.","python":sys.version,"elapsed_seconds":time.monotonic()-started,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n");print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":main()
