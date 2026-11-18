"""Certify a normalized mirror-symmetric constant-width counterexample to the R203 65-line host."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx
from certify_ah_tapered_minkowski_closure import geometry

POSE_LIMIT=F("0.00001901")


def exact(x):
    q=F(x)
    return arb(q.numerator)/q.denominator


def decimal_hash(data):
    return sha256(data).hexdigest()


class Support:
    def __init__(self,spec,pi):
        assert F(spec["mirror_axis_degrees"])==30
        assert list(map(F,spec["curvature_piece_endpoints_degrees"]))==[0,20,F(45,2),F(115,2),60,90]
        dev=spec["curvature_deviation_from_half"]
        assert dev==["1/2","solve_zero_cos_moment","-1/2","-43811/200000","1/2"]
        self.breaks=list(map(F,spec["curvature_piece_endpoints_degrees"]))
        self.angle_unit=pi/180
        angles=[exact(w)*self.angle_unit for w in self.breaks]
        self.angles=angles
        s=[w.sin() for w in angles]
        c=[w.cos() for w in angles]
        q=[exact(dev[k]) if k!=1 else arb(0) for k in range(len(dev))]
        m=[b-a for a,b in zip(s,s[1:])]
        q[1]=-sum(p*v for p,v in zip(q,m,strict=True))/m[1]
        assert q[1]>exact(F(-1,2)) and q[1]<exact(F(1,2))
        assert sum(p*v for p,v in zip(q,m,strict=True)).contains(0)
        self.dev=q

        # On a piece f=h-1/2 is A cos(q)+B sin(q)+g.
        ac=arb(0)
        bc=arb(0)
        self.rows=[]
        for j,v in enumerate(q):
            self.rows.append((ac-v*c[j],bc-v*s[j],v))
            ac+=v*(c[j+1]-c[j])
            bc+=v*(s[j+1]-s[j])
        mid=30*self.angle_unit
        k=self.select(F(30))
        x,y,z=self.rows[k]
        self.centre=-(x*mid.cos()+y*mid.sin()+z)/mid.cos()
        self.rows=[(a+self.centre,b,c) for a,b,c in self.rows]
        assert self.at_exact(F(30)).contains(0) and self.at_exact(F(90)).contains(0)

    def select(self,m):
        for j,(a,b) in enumerate(zip(self.breaks,self.breaks[1:])):
            if a<=m<=b:
                return j
        raise ValueError(m)

    def at_exact(self,q):
        a,b,c=self.rows[self.select(q)]
        s,t=(exact(q)*self.angle_unit).sin_cos()
        return a*t+b*s+c

    def along(self,offset,sign,midpoint):
        a,b,c=self.rows[self.select(offset+sign*midpoint)]
        s,t=(exact(offset)*self.angle_unit).sin_cos()
        return a*t+b*s,sign*(b*t-a*s),c

    def compatibility_row(self,left,right):
        # Mirror about 30 degrees reduces [0,60] to [0,30].
        # On this half, F=f(30-phi)-f(30+phi)-f(90-phi).
        m=(left+right)/2
        a=self.along(F(30),-1,m)
        b=self.along(F(30),1,m)
        c=self.along(F(90),-1,m)
        return tuple(x-y-z for x,y,z in zip(a,b,c,strict=True))

    def folded_angle(self,u):
        x,y=u
        q=arb.atan2(y,x)/self.angle_unit-30
        while q>180:q-=360
        while q< -180:q+=360
        q=abs(q)
        sign=1
        if q>90:
            q=180-q
            sign=-1
        assert q>0 and q<90
        return q,sign

    def support_direction(self,u):
        q,sign=self.folded_angle(u)
        s,c=(q*self.angle_unit).sin_cos()
        vals=[]
        for j,(l,r) in enumerate(zip(self.breaks,self.breaks[1:])):
            if not (q<exact(l) or q>exact(r)):
                a,b,h=self.rows[j]
                vals.append(a*c+b*s+h)
        assert vals
        out=vals[0]
        for v in vals[1:]:out=out.union(v)
        return exact(F(1,2))+sign*out


class Polygon:
    """Match the two trigonometric moments in each nonbinary-curvature cell."""
    def __init__(self,base):
        self.base=base
        hi=exact(F(1,2))
        lo=-hi
        self.breaks=[base.angles[0]]
        self.dev=[]
        self.error=arb(0)
        self.replacements=[]
        for j,(a,b) in enumerate(zip(base.angles,base.angles[1:])):
            if j not in (1,3):
                self.dev.append(base.dev[j])
                self.breaks.append(b)
                continue
            period=b-a
            rho=base.dev[j]+hi
            mid=(a+b)/2
            half=(rho*(period/2).sin()).asin()
            left,right=mid-half,mid+half
            assert a<left<right<b
            self.breaks.extend((left,right,b))
            self.dev.extend((lo,hi,lo))
            # The centred on-arc matches rho times each of the two moments.
            for fn in ("sin","cos"):
                same=(getattr(right,fn)()-getattr(left,fn)())-rho*(getattr(b,fn)()-getattr(a,fn)())
                assert same.contains(0)
            self.error=max(self.error,(period**2/2).upper())
            self.replacements.append(tuple(map(str,(left,right))))
        self.rows=[]
        x=base.centre
        y=arb(0)
        for l,r,q in zip(self.breaks[:-1],self.breaks[1:],self.dev,strict=True):
            s,c=l.sin_cos()
            t,d=r.sin_cos()
            self.rows.append((x-q*c,y-q*s,q))
            x+=q*(d-c)
            y+=q*(t-s)
        assert len(self.rows)==9

    def support_direction(self,u):
        q,sign=self.base.folded_angle(u)
        t=q*self.base.angle_unit
        s,c=t.sin_cos()
        values=[a*c+b*s+d for l,r,(a,b,d) in zip(self.breaks[:-1],self.breaks[1:],self.rows,strict=True)
                if not (t<l or t>r)]
        assert values
        out=values[0]
        for v in values[1:]:
            out=out.union(v)
        return exact(F(1,2))+sign*out


def compatible_angles(support):
    points=[F(0),F("7.5"),F(10),F("27.5"),F(30)]
    res=[]
    for index,(l,r) in enumerate(zip(points,points[1:])):
        a,b,c=support.compatibility_row(l,r)
        best=None
        for k in range(64):
            p=l+(r-l)*F(k,64);q=l+(r-l)*F(k+1,64)
            x=(exact(p)*support.angle_unit).union(exact(q)*support.angle_unit)
            si,co=x.sin_cos()
            # At zero the exact normalization gives F(0)=0.
            val=-a*si+b*co if index==0 else a*co+b*si+c
            assert val<0,(l,r,k,str(val))
            best=val.upper() if best is None else max(val.upper(),best)
        if index>0:
            assert best<exact(F("-.0189"))
        res.append({"domain_degrees":[str(l),str(r)],"quantity":"derivative per radian" if index==0 else "F",
                    "upper":str(best),"closed_cells":64})
    return res


def parallel_cut(env,shift):
    g=env["g"]
    u,h=env["u0"],env["h0"]+exact(shift)
    slack=[h-g.dot(u,v) for v in env["xx"]]
    retained=next(j for j,v in enumerate(slack) if v>0)
    assert all(v<0 for v in slack[:retained]) and all(v>0 for v in slack[retained:])
    assert h-g.dot(u,env["ns"]["z"])<0 and retained<=env["keep"]
    _,normal,offset=env["rows"][retained-1]
    edge=g.line_pair(u,normal,h,offset)
    assert g.dot(g.minus(edge,env["xx"][retained-1]),g.minus(env["xx"][retained],edge))>0
    ns=env["ns"]
    old_parameter=exact(F(".95"))
    param=old_parameter+exact(shift)/g.dot(u,ns["D"])
    assert 0<param<old_parameter
    point=g.plus(ns["z"],g.scale(param,ns["D"]))
    between=env["xx"][retained:env["keep"]]
    vertices=[point,edge,*between,env["crosspoint"],env["W"]]
    twice=sum(g.cross(x,y) for x,y in zip(vertices,vertices[1:]+vertices[:1],strict=True))
    area=abs(twice)/2
    assert area>0
    return area,retained


def repair_rows(env,cloud):
    out=[]
    rows={label:(n,b) for label,n,b in env["active"]}
    # These initially inactive rows can become active after a parallel shift.
    rows.update({label:(n,b) for label,n,b in env["rows"]})
    for label,(n,b) in rows.items():
        if label in env["hexnames"] or label.startswith(env["sprague"]) or label in ("taper_64","deep"):
            continue
        u=env["preimage"](n,240,1)
        clearance=b-cloud.support_direction(u)
        assert clearance>exact(F(".00001901")),(label,str(clearance))
        out.append((label,clearance.lower()))
    return out


def main():
    p=ArgumentParser(description=__doc__)
    p.add_argument("computation",type=Path)
    p.add_argument("input",type=Path)
    args=p.parse_args()
    env,digest=geometry(args.computation)
    ctx.prec=192
    pi=arb.pi()
    raw=args.input.read_bytes()
    spec=json.loads(raw)
    host={row[0]:row[1:] for row in env["active"]}
    data=Support(spec,pi)
    whole=compatible_angles(data)
    polygon=Polygon(data)
    # The endpoint chamber is undisturbed: q=30-phi,30+phi,90-phi
    # all lie outside the interiors of the two replaced intervals.
    assert 3*polygon.error<exact(F(".003"))
    g=env["g"]
    for label,i in (("I-E3",3),("K-C3",5)):
        u,b=host[label]
        assert all((a-z).contains(0) for a,z in zip(u,env["tilt"][i]))
        assert (b-exact(F(1,2))).contains(0)
    rows=[]
    polygon_rows=[]
    assert [k for k,_ in spec["rotation_faces"]]==list(range(0,360,60))
    for k,label in spec["rotation_faces"]:
        assert label in ("I-E3","K-C3","deep")
        n,b=host[label]
        u=env["preimage"](n,k,1)
        assert g.norm(u).contains(1)
        violation=data.support_direction(u)-b
        assert violation>exact(POSE_LIMIT),(k,label,str(violation))
        rows.append({"rotation_degrees":k,"true_face":label,"violation_lower":str(violation.lower()),
                     "violation_display":float(violation.lower())})
        pviolation=polygon.support_direction(u)-b
        assert pviolation>exact(POSE_LIMIT),(k,label,str(pviolation))
        polygon_rows.append({"rotation_degrees":k,"true_face":label,"violation_lower":str(pviolation.lower()),
                             "violation_display":float(pviolation.lower())})
    # An explicit support-value difference rules out the previously studied 3-fold class.
    asym=data.at_exact(F(0))-data.at_exact(F(60))*(-1)
    assert abs(asym)>exact(F(".001"))
    # Determine the original mirror-body repair cost in the only orientation
    # permitted by both unchanged slanted faces.
    other=repair_rows(env,data)
    other_polygon=repair_rows(env,polygon)
    edge,margin=min(other,key=lambda row:row[1])
    assert data.support_direction(env["preimage"](host["deep"][0],240,1))-host["deep"][1]<exact(F(".00001902"))
    pshift=F(".00002291")
    assert polygon.support_direction(env["preimage"](host["deep"][0],240,1))-host["deep"][1]<exact(pshift)
    a_lower,ear_index=parallel_cut(env,F(".0000229072"))
    a_upper,_=parallel_cut(env,pshift)
    target=env["oldarea"]-exact(F(".00001"))
    assert env["area"]+a_lower>target
    print(json.dumps({
        "complete":True,"fixed_host_65_line_counterexample":True,
        "domain":"all rigid isometries, with reflection reduced by the exact mirror",
        "input_sha256":decimal_hash(raw),"derived_host_geometry_sha256":digest,
        "mirror_axis_degrees":30,"curvature_deviation":[str(v) for v in data.dev],
        "centre_support_deviation":str(data.centre),
        "decomposition":"h=1/2+f; antiperiodic f, mirror-even around 30 degrees",
        "angular_side":"F(phi)<0 for every 0<phi<60 degrees",
        "compatibility_closed_cells_half_chamber":whole,
        "six_rotations":rows,"six_compatible_rotations_violation_lower":str(POSE_LIMIT),
        "support_difference_rotation_120":str(asym),
        "original_bounded_repair_at_240":{"sufficient_parallel_deep_AH_shift":"1902/100000000",
                                         "other_rows":len(other),"weakest_other_edge":edge,
                                         "strict_other_margin_lower":str(margin)},
        "parallel_repair_for_polygon":{
            "necessary_parallel_shift_lower":"28634/1250000000",
            "sufficient_parallel_shift":str(pshift),"other_rows":len(other_polygon),
            "first_retained_chain_vertex_at_lower":ear_index,
            "added_area_lower":str(a_lower.lower()),"added_area_at_sufficient_shift":str(a_upper),
            "B_grade_area_deficit_lower":str((env["area"]+a_lower-target).lower())},
        "binary_curvature_Reuleaux_polygon":{
            "full_period_unit_arc_count":13,"same_compatibility_roots":True,
            "Hausdorff_error_upper":str(polygon.error),
            "changed_q_intervals_degrees":[["20","45/2"],["115/2","60"]],
            "centred_unit_radius_intervals_radians":polygon.replacements,
            "six_rotations":polygon_rows},
    },indent=2))


if __name__=="__main__":
    main()
