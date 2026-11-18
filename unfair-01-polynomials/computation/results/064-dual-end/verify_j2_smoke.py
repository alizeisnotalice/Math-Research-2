"""Direct polynomial division checks for the C++ meet-in-the-middle program."""

import json
from pathlib import Path
import subprocess
import sys

import sympy as sp

Z=sp.Symbol("z")


def poly(values):
    return sp.Poly.from_list(values,Z,modulus=163)


def direct(degree,q,fixture):
    j=poly([1,6,121,106,68,115,117,q])
    modulus=j*j
    def image(mask,include_endpoints):
        c=[0]*(degree+1)
        if include_endpoints:
            c[:12]=[1]*7+[0]+[1]*4
            c[-1]=1
        for i in range(degree-12):
            c[i+12]=(mask>>i)&1
        return poly([value*(-1)**i for i,value in enumerate(c)]).rem(modulus)
    target=poly([0]) if fixture is None else image(fixture,False)
    return [mask for mask in range(1<<(degree-12)) if image(mask,fixture is None)==target]


def main():
    assert __debug__
    binary=str(Path(sys.argv[1]).resolve())
    cases=[(13,0,None),(16,52,None),(18,162,None),(20,0,37),
           (20,19,255),(18,5,0),(21,104,301)]
    report=[]
    for degree,q,fixture in cases:
        args=[binary,str(degree),str(q)]
        if fixture is not None:
            args.append(str(fixture))
        actual=json.loads(subprocess.check_output(args,text=True))
        expected=direct(degree,q,fixture)
        assert sorted(actual["matches"])==expected
        assert actual["left_patterns"]*actual["right_patterns"]==1<<(degree-12)
        if fixture is not None:
            assert fixture in actual["matches"]
        report.append(actual)
    print(json.dumps({"status":"verified","arithmetic_checks":report}))


if __name__=="__main__":
    main()
