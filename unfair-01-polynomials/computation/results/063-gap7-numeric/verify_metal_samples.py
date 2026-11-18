"""Check saved Swift Float states against an independent exact convolution."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import struct


def integer_coefficients(values):
    ratios=[struct.unpack("f",struct.pack("f",q))[0].as_integer_ratio() for q in values]
    scale=max(d for _,d in ratios)
    assert all(scale%d==0 for _,d in ratios)
    return [n*(scale//d) for n,d in ratios],scale


def main():
    assert __debug__
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    opts=parser.parse_args()
    raw=opts.input.read_bytes()
    data=json.loads(raw)
    assert data["status"]=="finished_numeric_search" and data["completed"]==data["steps"]
    assert data["device"] and len(data["cases"])==16
    m,n=data["m"],data["n"]
    output=[]
    for row in data["cases"]:
        a,da=integer_coefficients(row["A"])
        b,db=integer_coefficients(row["B"])
        assert (len(a),len(b))==(m+1,n+1)
        assert a[0]==a[-1]==da and b[0]==b[-1]==db
        assert all(0<=x<=da for x in a) and all(0<=x<=db for x in b)
        assert all(a[i]==0 for j in (m,n) for gap in (1,2,5) for i in [j-gap] if i<=m)
        assert all(b[i]==0 for j in (m,n) for gap in (1,3,4) for i in [j-gap] if i<=n)
        values=[0]*(m+n+1)
        for i,u in enumerate(a):
            for j,v in enumerate(b):
                values[i+j]+=u*v
        scale=da*db
        assert all(abs(Fraction(values[k],scale)-1)<Fraction(1,10**6) for k in range(7))
        assert values[7]<scale//10**6
        targets=[int(2*v>=scale or k in (8,9,10,11,m,n)) for k,v in enumerate(values)]
        error=max(abs(Fraction(values[k],scale)-targets[k]) for k in range(8,len(values)))
        assert abs(float(error)-row["residual"])<2e-5
        output.append(str(error))
    report={"status":"verified_rounded_samples","samples":len(output),"sha256":sha256(raw).hexdigest(),
            "device":data["device"],"residuals":output}
    with opts.out.open("x") as f:
        json.dump(report,f,indent=2)
        f.write("\n")
    print(json.dumps({"status":report["status"],"samples":len(output),"device":data["device"]}))


if __name__=="__main__":
    main()
