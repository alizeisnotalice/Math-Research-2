"""Bound every near-minimal five-body convex host in Hausdorff distance from the unique minimizing hull."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,ctx

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("contact",type=Path)
    ap.add_argument("--cells",type=int,default=12)
    a=ap.parse_args();source=a.contact.resolve();sys.path.insert(0,str(source))
    import certify_five_body_contact_basin as b
    from four_contact_taylor import powers
    ctx.prec=192
    basin=source/"five-body-contact-basin.json";record=json.loads(basin.read_bytes())
    assert sha256(basin.read_bytes()).hexdigest()=="282fafa2a4d37f72871c082b9e3e3e71c19fdb9ee4889cfa1a08fef19dcd150a"
    assert F(record["normal_halfwidth"])==F(3,1000000)and record["degree"]==9
    radius=b.q(2*F(record["normal_halfwidth"]));cx,cy=map(b.q,record["center"])
    b.Model.degree=record["degree"]
    u=b.Model({(0,0):cx,(1,0):radius});v=b.Model({(0,0):cy,(0,1):radius});pose=b.model(u,v)[1]
    def jac(z,axis,x,y):
        xp,yp=powers(x,z.degree),powers(y,z.degree)
        val=arb(0)
        for (m,n),c in z.poly.items():
            if axis==0 and m:val+=m*c*xp[m-1]*yp[n]
            if axis==1 and n:val+=n*c*xp[m]*yp[n-1]
        margin=1-max(x.abs_upper(),y.abs_upper());assert margin>0
        return arb(val,(z.error/margin).abs_upper())/radius
    def matrix_norm(j):
        aa=sum(v[0]**2 for v in j);cc=sum(v[1]**2 for v in j);bb=sum(v[0]*v[1]for v in j)
        radicand=(aa-cc)**2+4*bb**2
        return ((aa+cc+radicand.abs_upper().sqrt())/2).abs_upper().sqrt()
    count=a.cells;assert count>0
    rows={n:[arb(0),arb(0),arb(0)]for n in (3,5,7,9)}
    for i in range(count):
        for j in range(count):
            unit=lambda k:arb(b.q(F(2*k+1-count,2*count)),b.q(F(1,2*count)))
            x=unit(i);y=unit(j);J=[[jac(z,d,x,y)for d in (0,1)]for z in pose]
            for body,n in enumerate((3,5,7,9)):
                k=3*body
                tlip=matrix_norm(J[k:k+2])
                rot=arb(0)if body==0 else (J[k-1][0]**2+J[k-1][1]**2).abs_upper().sqrt()
                ell=(b.pi/(2*n)).tan()/2 if body else arb(0)
                for idx,value in enumerate((tlip,rot,tlip+ell*rot)):
                    rows[n][idx]=max(rows[n][idx],value.abs_upper())
    L=max(row[2]for row in rows.values())
    stat=source/"five-body-contact-stationary.json";s=json.loads(stat.read_bytes())
    assert sha256(stat.read_bytes()).hexdigest()=="1ad6912c1939657959fda04e6fea3312688516343ddadde9803ce70f522c377a"
    U=b.q(".83729500622102525");M=arb(s["area_ball"]);gap=(U-M).abs_upper()
    R=b.q(".703");r=b.q(".5")
    point_area=b.pi*r*r+r*(R*R-r*r).sqrt()-r*r*(r/R).acos()
    assert point_area>U
    host_sq=(4*R*R-1).sqrt()
    alpha=b.q(".05");beta=b.q(".0000167")
    C=(host_sq+L*L/alpha).abs_upper().sqrt()
    # Root perimeter deficit alone gives C sqrt(e) throughout this energy range.
    assert 2*gap.sqrt()<C*beta
    raw=C*gap.sqrt()+gap/beta;single=C*gap.sqrt()
    C0=b.q("5.381635");D0=b.q("59880.24")
    assert C<C0 and 1/beta<D0 and C0*gap.sqrt()<b.q(".00000001376648")
    assert C0*gap.sqrt()+D0*gap<b.q(".00000001376687")
    print(json.dumps(dict(rounded_sqrt="5.381635",rounded_linear="59880.24",rounded_threshold_bound="0.00000001376648",target=str(F(".83729500622102525")),threshold_decimal=".83729500622102525",
        denominator_normal="0.05",denominator_transverse="0.0000167",cells=count**2,host_circumradius="0.703",
        hull_of_disk_and_radius_point_area=str(point_area),per_body={n:dict(translation_op_upper=str(q[0]),angle_gradient_norm_upper=str(q[1]),
            Hausdorff_per_unit_normal_upper=str(q[2]))for n,q in rows.items()},chart_upper=str(L),C_sqrt=str(C),C_linear=str(1/beta),
        threshold_area_excess_upper=str(gap),threshold_raw_Hausdorff_upper=str(raw),threshold_Hausdorff_upper=str(single),
        old_basin_sha256=sha256(basin.read_bytes()).hexdigest(),stationary_sha256=sha256(stat.read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
if __name__=="__main__":main()
