"""Check an outer model for every width-one body with a uniform curvature-radius floor."""
from hashlib import sha256
from itertools import pairwise
from pathlib import Path
import json
import zipfile

from flint import arb, ctx


def scalar_product(p, q):
    return p[0] * q[0] + p[1] * q[1]


def verify():
    ctx.prec = 224
    root = Path(__file__).resolve().parent
    with zipfile.ZipFile(root / "degree-nine.zip") as archive:
        raw = archive.read("host.json")
    digest = sha256(raw).hexdigest()
    assert digest == "a69bf398326ebc2e785cdecbcf0386853dd8b0622eea5d6ac06e09b259f0194f"
    host = json.loads(raw)
    denominator = host["denominator"]
    vertices = host["vertices"]

    degrees = arb.pi() / 180
    sigma = arb(1549894) / 10**6
    phi = arb(106) / 10**6
    angles = [0, 60, 120, 180, 210 + sigma, 240, 300, 330 + sigma]
    # Rotating the entire base polygon rotates each normal by phi.
    normals = []
    for angle in angles:
        sine, cosine = ((angle + phi) * degrees).sin_cos()
        normals.append((cosine, sine))
    outer = []
    for a, b in pairwise(normals + normals[:1]):
        divisor = 2 * (1 + scalar_product(a, b))
        assert divisor > 0
        outer.append(((a[0] + b[0]) / divisor, (a[1] + b[1]) / divisor))

    t = arb(441834) / 10**6
    shift = (arb(36719) / 10**12, -arb(68456) / 10**12)
    worst = None
    for a, b in pairwise(vertices + vertices[:1]):
        u = (b[1] - a[1], a[0] - b[0])
        length = arb(scalar_product(u, u)).sqrt()
        assert length > 0
        height = arb(scalar_product(a, u)) / denominator
        for vertex in outer:
            support = scalar_product(u, shift) + (1 - t) * length / 2 + t * scalar_product(u, vertex)
            slack = (height - support) / length
            assert slack > arb(8) / 10**12
            if worst is None or slack.lower() < worst:
                worst = slack.lower()
    return dict(host_sha256=digest, edge_count=len(vertices), base_vertex_count=len(outer),
                sigma_degrees="1.549894", rotation_degrees="0.000106",
                translation=["0.000000036719", "-0.000000068456"],
                t="0.441834", curvature_radius_interval=["0.279083", "0.720917"],
                minimum_slack_lower=str(worst))


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
