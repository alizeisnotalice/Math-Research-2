"""Transfer constants shared by the .83743 and .83744 signed screens."""
import json
from fractions import Fraction as F
from pathlib import Path

from certify_polygon_083731 import ar, arb, lower_radius, ORDERS
from certify_tomar_six_high import SourceTree
from certify_refined_material_contraction import SCALE, SHRINK, read


def verify(comp, tree_path):
    source=read(comp/"refined-material-contraction.json")
    f=F(source["point_factor"])
    assert F(source["curve_scale"])==SCALE and F(source["point_shrink"])==SHRINK
    s,r,R,h=SCALE,F(1,2),F(".703"),F(".01464")
    a=r-F(".069"); m=SHRINK*F(".528619")/(1-s)
    v=(1+s)*r*r+(1-s)*m*m
    derivative=1+s-(2*m*r+v)/(4*r*r)
    margin=4*(a-m)*s*s-(1+s)*h-max(0,(v-2*m*(r+h))/(2*r+h))
    cap=s*s*F(".837295")+s*r*h*h/R
    pi=arb.pi(); sr=ar(s*r); rr=ar(R)
    assert (pi-(sr/rr).acos())*sr**2+sr*(rr**2-sr**2).sqrt()>ar(cap)
    assert R*R<2*(s*r)**2 and m<a<r and derivative>0 and margin>0 and cap>F(".83744577")
    loss=F(source["exceptional_loss"])
    tree=SourceTree(tree_path)
    rows=[]
    for code in ("083743","083744"):
        seed=read(comp/f"tomar-{code}-seed.json")
        G=F(seed["target_integer"],2**57)
        assert F(seed["target"])<G<cap
        assert all(F(x["integer_radius"],2**44)<F(".069") and F(x["integer_area_lower"],2**57)>G for x in seed["radii"])
        for j,n in enumerate(ORDERS[1:]):
            k=3*j
            d=min(F(tree.root["hi"][k]),F(tree.root["hi"][k+1]),-F(tree.root["lo"][k]),-F(tree.root["lo"][k+1]))
            ins=ar(s)*lower_radius(n,pi)
            assert pi*(sr**2+ins**2)/2+ar(d)*(sr+ins)>ar(G)
        high=F(seed["target"])+F(".000005")
        assert high-loss>G
        rows.append(dict(target=seed["target"],integer_goal=str(G),exceptional_required=str(high),
                         conditional_floor=str(f*f*G),display=float(f*f*G)))
    return dict(cap_height=str(h),common_cap_limit=str(cap),comparison_slack=str(margin),
                derivative=str(derivative),exceptional_loss=str(loss),conditional_rows=rows)


if __name__=="__main__":
    import sys
    print(json.dumps(verify(Path(__file__).resolve().parent,Path(sys.argv[1])),indent=2))
