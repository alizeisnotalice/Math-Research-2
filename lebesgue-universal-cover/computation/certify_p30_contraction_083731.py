"""Check remaining cap margin and diameter renormalization at R136."""
from fractions import Fraction as F
from hashlib import sha256
import json

from certify_finite_point_leaf import U
from certify_nine import ar, arb, ctx
from certify_polygon_083731 import AREA, GOAL, HERE, LAMBDA, ORDERS, lower_radius


def main():
    ctx.prec=160
    source=(HERE/"finite-point-leaf-input.json").read_bytes()
    spec=json.loads(source)
    denom=spec["denominator"]
    points=[(F(x,denom), F(y,denom)) for x,y in spec["points"]]
    p30_diameter_sq=max((x-u)**2+(y-v)**2 for j,(x,y) in enumerate(points) for u,v in points[:j])
    assert p30_diameter_sq < 1
    center_x=F("-.01007")
    point_radius=F(".528619")
    point_radius_sq=max((x-center_x)**2+y*y for x,y in points)
    assert point_radius_sq < point_radius**2
    rho=F(713,10**9)
    delta=rho*point_radius
    s=LAMBDA
    eps=1-s
    r=F(1,2)
    a=F(".42264973")
    pi=arb.pi()
    assert all(lower_radius(n,pi) > ar(a) for n in ORDERS)
    c=delta/eps
    assert c<a<r

    R=F(".7022")
    R_a=ar(R)
    r_a=ar(s*r)
    spike=(pi-(r_a/R_a).acos())*r_a**2+r_a*(R_a**2-r_a**2).sqrt()
    assert spike > ar(F(".8373185106")) > ar(F(GOAL,AREA))
    # Winning cap direction: x_old > b+d, and x_old <= R+delta.
    C=2**20
    b_round=(R+delta)*F(4,U-2)+F(1,C)
    assert b_round < F(".000002")
    inflated_x=F(".702203")
    assert R+delta+b_round < inflated_x
    gap=F(1,4000)+F(1,C**2)
    depth=F(".01571")
    d0=F(int(depth*C),C)
    assert d0**4 > gap**2*(4*inflated_x**2-1)

    # Compare the old hat at e to the new one with both point and wall shifted.
    polynomial=(1+s)*r**2+c*c*eps
    derivative_lower=(1+s)-(2*c*r+polynomial)/(4*r*r)
    assert derivative_lower>0
    upper_factor=(1+s)*depth+max(F(0),(polynomial-2*c*(r+depth))/(2*r+depth))
    lower_factor=4*(a-c)*s*s
    assert lower_factor-upper_factor>F(".00122")

    reduced_diameter_sq=max(s*s,(1-rho)**2*p30_diameter_sq)
    assert reduced_diameter_sq==(1-rho)**2*p30_diameter_sq
    floor=F(GOAL,AREA)/reduced_diameter_sq
    rational_floor=F(GOAL,AREA)/(1-rho)**2
    assert floor>F(".83731119527") and rational_floor>F(".837311194")
    out={
        "P30_input_sha256":sha256(source).hexdigest(),
        "integer_area_floor":str(F(GOAL,AREA)),
        "core_scale":str(s), "P30_contraction":str(rho), "P30_center_x":str(center_x),
        "P30_radius_squared_exact":str(point_radius_sq), "P30_radius_upper":str(point_radius),
        "P30_diameter_squared_exact":str(p30_diameter_sq),
        "point_projection_loss_upper":str(delta), "min_centered_support_lower":str(a),
        "radius_clip":str(R), "radius_clip_area_lower":str(spike),
        "unit_support_rounding_upper":str(b_round), "total_shifted_projection_upper":str(inflated_x),
        "old_gap_upper":str(gap), "old_winning_cap_depth_upper":str(depth), "dyadic_test_depth":str(d0),
        "comparison_margin":str(lower_factor-upper_factor),
        "reduced_common_diameter_squared":str(reduced_diameter_sq), "normalized_area_floor":str(floor),
        "rational_normalization_factor":str(1/(1-rho)),
        "normalized_core_scale_rational":str(s/(1-rho)),
        "rational_normalization_area_floor":str(rational_floor)
    }
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
