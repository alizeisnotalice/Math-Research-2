"""Area floor using the infeasible fixed-EH endpoint and the compulsory AH chord."""
from fractions import Fraction as F
import json
from pathlib import Path
import audit_optimized_eh_chord as geo
from certify_eh_joint_parameters import specialize, gain_endpoint

here = Path(__file__).resolve().parent
def read(name):
    return json.loads((here/name).read_bytes())

sigma, endpoint = F("1.5497"), "0.00055222"
specialize(sigma)
gain = geo.arb(geo.trial(endpoint,tight=True)["gain_area"])
obstructions = read("eh-no-interior-roots-552p22-certificate.json")
assert any(F(row["delta"])==F(endpoint) and F(row["sigma"])==sigma
           and row["only_d6_poses_in_P"]==["identity","Q2"] for row in obstructions)
old = read("gibbs-1p5497-refined-18-128.json")
assert F(old["sigma_degrees"])==sigma
tail = read("ah-compulsory-chord-area.json")
floor = F(old["area_interval"][0])-gain_endpoint(gain.upper())-F(tail["maximum_AH_cap_area_interval_exact"][1])
print(json.dumps(dict(EH_limit=endpoint,EH_gain_ceiling=str(gain_endpoint(gain.upper())),
                      class_area_lower=str(floor),down_decimal="0.%015d"%(floor*10**15//1)),indent=2))
