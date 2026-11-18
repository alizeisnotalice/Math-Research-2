"""Size the five uniform samplings for an all-rational six-input certificate."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

from certify_nine import S, append_point, ar, arb, arc_body, ctx

HERE = Path(__file__).resolve().parent
PI = F(22,7)
T = F('.837295006222')


def area(e):
    return (T - PI*e**2)/(1 + 2*e/(F(1,2)-e))


def error(n, m):
    return PI**2/(8*(n*m)**2) + F(12,S) if n else PI**2/(4*m*m) + F(12,S)


def size(n, target):
    def enough(m):
        e = error(n, m)
        return e < F(1,2) and area(e) > target
    a, b = 0, 8
    while not enough(b):
        b *= 2
    while b-a > 1:
        mid = (a+b)//2
        if enough(mid):
            b = mid
        else:
            a = mid
    return b


def main(out):
    out.mkdir(exist_ok=True)
    ctx.prec = 160
    pi = arb.pi()
    assert pi < ar(PI)
    fixed = json.loads((HERE/'finite-point-leaf-input.json').read_bytes())
    assert sha256((HERE/'finite-point-leaf-input.json').read_bytes()).hexdigest() == 'b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc'
    results = []
    for label in ('.83729', '.837295'):
        target = F(label)
        groups = []
        errors = []
        for n in (0,3,5,7,9):
            m = size(n, target)
            points = []
            if n:
                arc_body(points,m,pi,n,arb(0),arb(0),arb(0))
            else:
                for j in range(m):
                    sn, cs = (2*pi*j/m).sin_cos()
                    append_point(points,cs/2,sn/2)
            points = sorted(set(points))
            assert len(points) == max(n,1)*m
            groups.append(dict(n=n, denominator=S+10, points=points))
            errors.append(error(n,m))
        groups.append(dict(n='P30',denominator=fixed['denominator'],points=fixed['points']))
        path = out/('six-rational-'+label[1:]+'.json')
        raw = (json.dumps(groups,separators=(',',':'))+'\n').encode()
        path.write_bytes(raw)
        counts = [len(g['points']) for g in groups]
        results.append(dict(target=label,counts=counts,total=sum(counts),error=str(max(errors)),
                            area_lower=str(area(max(errors))),sha256=sha256(raw).hexdigest()))
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
