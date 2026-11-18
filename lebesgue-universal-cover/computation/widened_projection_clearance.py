"""Transfer the exact projected scale bounds to a Euclidean outer parallel body."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

root = Path(__file__).resolve().parent
cert = json.loads((root/"widened-projected-obstruction-certificate.json").read_bytes())
host = root/"r112-widened-host-certificate.json"
assert sha256(host.read_bytes()).hexdigest() == cert["source_sha256"][host.name]
planes = json.loads(host.read_bytes())["rational_planes"]
assert all(c>0 and 4*c*c >= x*x+y*y for _,x,y,c in planes)
distance = F(211,10**7)
scale = 1/(1+2*distance)
squared = [F(row["max_squared_radius_exact"]) for row in cert["branches"]]
assert len(squared)==2 and max(squared)<scale**2
print(json.dumps(dict(centered_inradius_lower="1/2", clearance=str(distance),
                      required_squared_scale=str(scale**2), both_chiralities_excluded=True),indent=2))
