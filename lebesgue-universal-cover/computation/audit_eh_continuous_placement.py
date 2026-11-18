#!/usr/bin/env python3
"""Construct an inscribed polygon for Gibbs's chord cut and check rigid placements."""
import argparse
import json
from pathlib import Path

from audit_optimized_eh_chord import (
    HALF, ONE, arb, circle_line, circle_pair, cross, dot, line_pair, minus, norm,
    plus, rot, scale, setup, unit,
)


def project(point, direction):
    return minus(point,scale(dot(point,direction)-HALF,direction))


def polygon(delta, subdivisions=1):
    n,m,c3,b3,_,u,_,te,_=setup()
    b1,f1,d1=(line_pair(n[i],n[j]) for i,j in ((0,1),(2,3),(4,5)))
    e3,f3=line_pair(n[4],m[3]),line_pair(n[3],m[2])
    g=circle_line(f3,n[5],1)
    q,p,ip,h,k=(project(x,y) for x,y in ((e3,n[1]),(g,n[2]),(b3,m[3]),(c3,n[3]),(f3,m[5])))

    r3=arb(3).sqrt()
    c=(-e3[0],e3[1])
    midpoint=(-HALF/2,-r3/4)
    l0=(r3/2-HALF,(r3/2-HALF-1)/r3)
    z=circle_pair(l0,midpoint,-1)
    y=(arb(0),c[1]+(1-c[0]**2).sqrt())
    v=minus(z,c)
    v2=dot(v,v)
    ta=plus(c,scale(1/v2,minus(v,scale((v2-1).sqrt(),rot(v)))))

    tangent=minus(ta,z)

    def paired(parameter):
        base=plus(z,scale(parameter,tangent))

        def inverse(normal,along):
            d=dot(normal,base)-HALF
            dd=dot(normal,tangent)
            root=(1-d*d).sqrt()
            value=minus(minus(base,scale(d,normal)),scale(root,along))
            derivative=plus(minus(tangent,scale(dd,normal)),scale(d*dd/root,along))
            return value,derivative

        l,dl=inverse((HALF,-r3/2),(r3/2,HALF))
        nn,dnn=inverse((-HALF,-r3/2),(-r3/2,HALF))
        v0=l[0]+HALF
        root=(1-v0*v0).sqrt()
        lx=(v0+r3*root)/2-HALF
        dlx=(1-r3*v0/root)*dl[0]/2
        a,b=(lx,(lx-1)/r3),minus(scale(2,midpoint),nn)
        da,db=(dlx,dlx/r3),scale(-1,dnn)
        step,dstep=minus(b,a),minus(db,da)
        size=dot(step,step)
        factor=(1/size-HALF/2).sqrt()
        dfactor=-dot(step,dstep)/(size*size*factor)
        value=minus(scale(HALF,plus(a,b)),scale(factor,rot(step)))
        derivative=minus(scale(HALF,plus(da,db)),plus(scale(dfactor,rot(step)),scale(factor,rot(dstep))))
        return value,derivative

    tp,_=paired(ONE)

    co,si=unit(arb(delta))
    radial=minus(u,c3)
    u1=plus(c3,minus(scale(co,radial),scale(si,rot(radial))))
    assert cross(minus(h,c3),minus(u1,c3))>0
    assert cross(minus(u1,c3),minus(u,c3))>0
    points={"B1":b1,"Q":q,"Y":y,"Ta":ta,"Z":z,"Tp":tp,"P":p,"F1":f1,"H":h,
            "U1":u1,"Te":te,"I":ip,"E3":e3,"D1":d1,"G":g,"K":k,"C3":c3}
    chord=minus(te,u1)
    assert all(cross(chord,minus(point,u1))>0 for name,point in points.items() if name not in ("U1","Te"))
    arc_centers={"Q":e3,"Y":c,"Tp":g,"H":c3,"Te":b3,"G":f3}
    boundary=list(points.items())
    edges=[]
    outer_edges=[]
    for (name,point),(name2,point2) in zip(boundary,boundary[1:]+boundary[:1]):
        edge=minus(point2,point)
        length=norm(edge)
        assert length>0
        normal=scale(-1/length,rot(edge))
        if name not in arc_centers and name!="Z":
            segment=(name+"-"+name2,normal,dot(normal,point))
            edges.append(segment)
            outer_edges.append(segment)
            continue

        if name in arc_centers:
            center=arc_centers[name]
            a,b=minus(point,center),minus(point2,center)
            assert dot(a,b)>0 and cross(a,b)>0

            def point_and_normal(parameter):
                radial=plus(scale(1-parameter,a),scale(parameter,b))
                normal=scale(1/norm(radial),radial)
                return plus(center,normal),normal
        else:
            def point_and_normal(parameter):
                value,derivative=paired(parameter)
                assert derivative[0]<0 and derivative[1]<0
                return value,scale(-1/norm(derivative),rot(derivative))

        outer_count=max(4,subdivisions)
        for i in range(outer_count+1):
            sample,normal=point_and_normal(arb(i)/outer_count)
            outer_edges.append((f"tangent_{name}_{i}",normal,dot(normal,sample)))

        piece=[point]+[point_and_normal(arb(i)/subdivisions)[0] for i in range(1,subdivisions)]+[point2]
        for index,(a,b) in enumerate(zip(piece,piece[1:])):
            displacement=minus(b,a)
            normal=scale(-1/norm(displacement),rot(displacement))
            label=name+"-"+name2 if subdivisions==1 else f"chord_{name}_{index}"
            edges.append((label,normal,dot(normal,a)))
    return edges,outer_edges


def rigid_check(edges, points, angle, x, y):
    cosine,sine=unit(arb(angle))
    images={name:plus((arb(x),arb(y)),plus(scale(cosine,p),scale(sine,rot(p))))
            for name,p in points.items()}
    slacks=[(offset-dot(normal,p),label,name) for label,normal,offset in edges for name,p in images.items()]
    slack,label,name=min(slacks,key=lambda row: row[0].lower())
    return {"angle":angle,"translation":[x,y],"contained":all(s>0 for s,_,_ in slacks),
            "minimum_signed_distance_slack":str(slack),"closest_edge":label,"point":name}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("witness",type=Path)
    parser.add_argument("delta")
    parser.add_argument("angle")
    parser.add_argument("x")
    parser.add_argument("y")
    args=parser.parse_args()
    points={name:tuple(map(arb,p)) for name,p in json.loads(args.witness.read_text())["points"].items()}
    print(json.dumps(rigid_check(polygon(args.delta)[0],points,args.angle,args.x,args.y),indent=2))
