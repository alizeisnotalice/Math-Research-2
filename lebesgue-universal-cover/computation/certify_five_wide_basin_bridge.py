"""Validate a raw parameter box against the two-angle, nine-transverse chart certificate."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from flint import arb,ctx

def main():
    if not __debug__:raise SystemExit("Assertions required.")
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("source",type=Path)
    p.add_argument("--tube",type=Path,default=Path(__file__).with_name("five-wide-contact-tube.json"))
    p.add_argument("--triangle",default=".0000005");p.add_argument("--others",default=".000030")
    a=p.parse_args();sys.path.insert(0,str(a.source.resolve()))
    import certify_five_body_contact_basin as b
    ctx.prec=192;row=json.loads(a.tube.read_bytes())
    assert row["complete"] and row["center"]==list(b.CENTERS) and row["halfwidth"]=="3/500000" and row["transverse"]==".00007"
    src=Path(__file__).with_name("certify_five_wide_contact_tube.py")
    assert row["code_sha256"]==sha256(src.read_bytes()).hexdigest()
    assert sum(Fraction(num,int(n)**2)for n,num in row["accepted_cells"].items())==1
    C=list(map(b.q,row["center"]));R=b.q(row["scale"]);D=b.q(row["transverse"])
    X=[b.interval(*v)for v in b.BOX]
    assert all(abs(z-c)<R/2 for z,c in zip(X,C,strict=True))
    b.Model.degree=row["degree"]
    u=b.Model({(0,0):C[0],(1,0):R});v=b.Model({(0,0):C[1],(0,1):R})
    raw=[arb(b.q(v),b.q(a.triangle if i<2 else a.others))for i,v in enumerate(b.RAW)]
    A=b.unit_root(b.add(raw[:2],b.regular(3,0)),b.r,1);B=b.unit_root(b.add(raw[:2],b.regular(3,1)),b.r,-1)
    assert A[0]>0 and B[1]>0
    x=((A[1]/A[0]).atan()-C[0])/R;y=(b.pi/2-(B[0]/B[1]).atan()-C[1])/R
    assert abs(x)<b.q(".5") and abs(y)<b.q(".5")
    ref=b.evaluate(b.model(u,v)[1],x,y)
    tol=max(abs(raw[i]-ref[i]).abs_upper()for i in range(2,11));assert tol<D
    true=b.model(*X)[1]
    assert all(abs(z-b.q(v))<b.q(a.triangle if i<2 else a.others)for i,(z,v)in enumerate(zip(true,b.RAW)))
    print(json.dumps(dict(tube_sha256=sha256(a.tube.read_bytes()).hexdigest(),raw_center=b.RAW,
        raw_halfwidth=[a.triangle]*2+[a.others]*9,normalized_angles=[str(x),str(y)],maximum_offset=str(tol)),indent=2))
if __name__=="__main__":main()
