#!/usr/bin/env python3
"""Outward polygon for the original K0. Arguments: computation directory, JSON output."""
import sys, json, math, hashlib, platform
from pathlib import Path
from fractions import Fraction as F
from flint import arb, ctx, __version__ as flint_version

def main(comp, target):
    ctx.prec=160
    sys.path.insert(0,str(comp))
    import audit_optimized_eh_chord as g
    from certify_ah_tapered_minkowski_closure import geometry
    old,_=geometry(comp)
    assert (g.SIGMA-arb('1.5497')*g.PI/180).contains(0)
    ns=old['ns']; n,m,c3,b3,_,u,_,te,j=g.setup(); A=arb(1)/2; lp=g.line_pair
    b1=lp(n[0],n[1]); f1=lp(n[2],n[3]); d1=lp(n[4],n[5]); e3=lp(n[4],m[3]); f3=lp(n[3],m[2]); G=g.circle_line(f3,n[5],1)
    def project(x,y): return g.plus(x,g.scale(A-g.dot(x,y),y))
    q=project(e3,n[1]); P=project(G,n[2]); I=project(b3,m[3]); H=project(c3,n[3]); K=project(f3,m[5]); c=(-e3[0],e3[1]); Y=(arb(0),c[1]+(1-c[0]**2).sqrt()); Ta=g.plus(ns['z'],ns['D']); Tp=old['q'](arb(1))
    # Branch selections for the algebraic endpoint identities in the proof.
    root3=arb(3).sqrt()
    midpoint=(-A/2,-root3/4)
    l0=(root3/2-A,(root3/2-A-1)/root3)
    along_l=(root3/2,A)
    along_n=(-root3/2,A)
    assert g.dot(g.minus(Ta,c),along_l)>arb('.1456152')
    assert g.dot(g.minus(ns['z'],l0),along_l)>arb('.1001027')
    assert g.dot(g.minus(ns['z'],midpoint),along_n)>arb('.3116240')
    assert abs(c[0]+A)<arb('.859772')
    assert g.dot(g.minus(Ta,c),g.neg(g.rot(ns['direct'])))>0
    assert all(offset-g.dot(normal,Ta)>0 for _,normal,offset in old['otherEH'])
    pts=dict(B1=b1,Q=q,Y=Y,Ta=Ta,Tp=Tp,P=P,F1=f1,H=H,J=j,I=I,E3=e3,D1=d1,G=G,K=K,C3=c3)
    seq=['B1','Q','Y','Ta','Tp','P','F1','H','J','I','E3','D1','G','K','C3','B1']; centers={'Q':e3,'Y':c,'Tp':G,'H':c3,'J':b3,'G':f3}
    # The two previously used sample subarcs must be strictly nested in the true arcs.
    co,si=g.unit(arb('0.0005522')); rad=g.minus(u,c3); u1=g.plus(c3,g.minus(g.scale(co,rad),g.scale(si,g.rot(rad))))
    ordering={}
    for label,center,a,b,z in [('H_U1_J',c3,H,u1,j),('J_Te_I',b3,j,te,I)]:
        x,y,w=(g.minus(p,center) for p in (a,b,z))
        v1,v2=g.cross(x,y),g.cross(y,w)
        assert v1>0 and v2>0 and g.dot(x,w)>0
        ordering[label]=[str(v1),str(v2)]
    pieces=[]; angle=arb(0); lines=[]; arcs={}; N=128
    for x,y in zip(seq,seq[1:]):
        a,b=pts[x],pts[y]
        if x in centers:
            c0=centers[x]; r=g.minus(a,c0); s=g.minus(b,c0); cr=g.cross(r,s); dt=g.dot(r,s)
            assert cr>0 and dt>0
            assert (g.dot(r,r)-1).contains(0) and (g.dot(s,s)-1).contains(0)
            angle+=arb.atan2(cr,dt); din,dout=g.rot(r),g.rot(s); arcs[x]=dict(cross=str(cr),dot=str(dt))
            for i in range(N+1):
                t=arb(i)/N; v=g.plus(g.scale(1-t,r),g.scale(t,s)); normal=g.scale(1/g.norm(v),v)
                lines.append((f'tangent_{x}_{i}',normal,g.dot(normal,c0)+1))
        else:
            din=dout=g.minus(b,a); normal=(din[1],-din[0]); lines.append((x+'-'+y,normal,g.dot(normal,a)))
        pieces.append((x,y,din,dout))
    joins={}
    for i,(x,y,din,dout) in enumerate(pieces):
        nxt=pieces[(i+1)%len(pieces)][2]; cr=g.cross(dout,nxt); dt=g.dot(dout,nxt); assert dt>0
        if y in {'Q','P','H','I','K'}: assert cr.contains(0)
        else: assert cr>0; angle+=arb.atan2(cr,dt)
        joins[y]=dict(cross=str(cr),dot=str(dt),tangent_identity=y in {'Q','P','H','I','K'})
    assert arb.pi()<angle<3*arb.pi() and angle.contains(2*arb.pi())
    assert len(old['otherEH'])==11
    lines.extend(old['otherEH'])
    assert not any(label.startswith('tangent_Z_') or label in {'U1-Te','Ta-Z'} for label,_,_ in lines)
    chord=next(h for label,_,h in lines if label=='Ta-Tp'); assert chord>0
    Q=10**12; V=10**12
    def quantize(label,u,h):
        x,y=[int(round(float(t.mid())*Q)) for t in u]
        hi=(h+abs(u[0]-arb(x)/Q)+abs(u[1]-arb(y)/Q)).upper()
        z=math.ceil(float(hi)*Q)+4
        while not arb(z)/Q>hi: z+=1
        return label,x,y,z
    rows=[quantize(*x) for x in lines]
    poly=[(F(-1),F(-1)),(F(1),F(-1)),(F(1),F(1)),(F(-1),F(1))]
    for _,a,b,c0 in rows:
        out=[]
        for j,(x,y) in enumerate(poly):
            px,py=poly[j-1]; aa=a*px+b*py-c0; bb=a*x+b*y-c0
            if (aa>0)!=(bb>0):
                t=aa/(aa-bb); out.append((px+t*(x-px),py+t*(y-py)))
            if bb<=0: out.append((x,y))
        poly=out; assert poly
    def half(seq):
        result=[]
        for p in seq:
            while len(result)>1 and (result[-1][0]-result[-2][0])*(p[1]-result[-1][1])-(result[-1][1]-result[-2][1])*(p[0]-result[-1][0])<=0: result.pop()
            result.append(p)
        return result
    def floor(x): return x.numerator//x.denominator
    points=[]
    for x,y in poly:
        ix=floor(x*V); iy=floor(y*V); points.extend((ix+dx,iy+dy) for dx in (0,1) for dy in (0,1))
    points=sorted(set(points)); hull=half(points)[:-1]+half(reversed(points))[:-1]
    a2=sum(x*v-y*u for (x,y),(u,v) in zip(hull,hull[1:]+hull[:1])); assert a2>0
    files=['audit_optimized_eh_chord.py','certify_ah_tapered_minkowski_closure.py','certify_ah_tapered_release.py','audit_eh_continuous_placement.py','certify_eh_joint_parameters.py','certify_eh_upper_chain.py','certify_ah_released_tangent.py']
    result={'scope':'outer polygon for K0 as defined in nonmonotone-ah-core-area.md; no universal-cover claim','sigma_degrees':'1.5497','precision_bits':ctx.prec,'subdivisions':N,'python':platform.python_version(),'flint':flint_version,'order_checks':ordering,'arcs':arcs,'joins':joins,'total_angle':str(angle),'excluded':['U1-Te','Ta-Z','tangent_Z_*'],'normal_grid':Q,'vertex_grid':V,'lines':len(rows),'exact_clip_vertices':len(poly),'outer_integer_vertices':hull,'double_area_integer':str(a2),'area_upper':str(F(a2,2*V*V)),'area_upper_decimal_display_only':str(float(F(a2,2*V*V))),'source_sha256':{f:hashlib.sha256((comp/f).read_bytes()).hexdigest() for f in files}}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('certified',len(rows),'halfplanes, area <=',result['area_upper'],'display',result['area_upper_decimal_display_only'])

if __name__=='__main__':
    if len(sys.argv)!=3: raise SystemExit('usage: certify_core.py COMPUTATION_DIR OUTPUT_JSON')
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
