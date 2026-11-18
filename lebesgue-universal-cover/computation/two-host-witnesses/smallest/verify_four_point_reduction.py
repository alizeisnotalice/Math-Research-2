"""Replay the two endpoint types with exact rationals and pin both original host definitions."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

PINS = {
    "certify_free_p30_upper.py": "459395f9445b58df6fce6108de73f567db1d3b33d5f518547fe493ad4129e118",
    "certify_nine.py": "610e94cf22e901453bc34f77afb8f0994d9e471432cbe20498fca3b1530c9541",
    "certify_seven_point_families.py": "b5c125b63b6215ecf177e97c908c181aaa76342c51bb31197668673f723562d9",
    "free-p30-reoptimized-pose.json": "21e6fceeb7f2033b69a6c2a7ae633abb26b6afd3621e32f6e6aa57f9b8770593",
    "seven-point-families-pose.json": "114b8e680bb0adaeb6453b1a87abe856616873e32dd6c01a946b3291e917ee8b",
}


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1]


def cross(p, q):
    return p[0] * q[1] - p[1] * q[0]


def check(name, vlen, offset, nfaces, ndisks):
    u = (F(1), F(0))
    v = (F(0), F(vlen))
    t = tuple(F(x) for x in offset)
    z = tuple(2 * x for x in t)
    alpha = cross(z, v) / cross(u, v)
    beta = cross(u, z) / cross(u, v)
    assert max(abs(alpha), abs(beta)) <= F(1, 2)
    faces = sum(2 * abs(x) == 1 for x in (alpha, beta))
    disks = 0
    for i, j in product((-1, 1), repeat=2):
        w = tuple((i * v[k] - j * u[k]) / 2 for k in range(2))
        assert dot(w, w) < 1
        far = tuple(z[k] + w[k] for k in range(2))
        assert dot(far, far) <= 1
        if dot(far, far) == 1:
            assert dot(t, far) > 0
            disks += 1
    assert (faces, disks) == (nfaces, ndisks)
    assert faces or disks
    return dict(
        name=name,
        u=list(map(str, u)),
        v=list(map(str, v)),
        t=list(map(str, t)),
        maximal_multiplier="2",
        parallelogram_equalities=faces,
        unit_side_equalities=disks,
    )


def main(comp, output):
    assert __debug__
    for name, digest in PINS.items():
        assert sha256((comp / name).read_bytes()).hexdigest() == digest
    p155 = json.loads((comp / "free-p30-reoptimized-pose.json").read_text())
    p175 = json.loads((comp / "seven-point-families-pose.json").read_text())["pose"]
    assert (len(p155), len(p175)) == (14, 18)
    assert p155[:2] == ["0.0081909227591130829", "-0.0088541638221488524"]
    assert p175[:2] == ["0.0032598755330460376", "-0.006855923199845153"]
    cases = [
        check("intersection reaches a diagonal endpoint", "1/5", ("0", "1/20"), 1, 0),
        check("one side reaches unit length", "4/5", ("3/20", "1/10"), 0, 1),
        check("two sides reach unit length", "10/13", ("11/52", "0"), 0, 2),
    ]
    result = dict(
        python=platform.python_version(),
        host_source_sha256=PINS,
        disk_radius="1/2",
        disk_center=["0", "0"],
        B3_angle="0",
        B3_shifts={"H155": p155[:2], "H175": p175[:2]},
        endpoint_examples=cases,
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("computation", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.computation, args.output)
