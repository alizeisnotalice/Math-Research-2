"""Check an explicit rational finite sampling using the six-body proof's integer slack."""
from fractions import Fraction as F
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
source = json.loads((HERE/'tomar-outer-083730-audit.json').read_bytes())
T = F('0.83730')
unit = 2**57
minimum = F(source['int_target'], unit)
assert minimum - T == F(source['smallest_score_accept_margin']) > 0
assert T + F(source['minimum_radial_margin']) > minimum
assert minimum < F('.838')
den = 10**22
n_disk, k_arc, n_smooth = 2**31, 2**31, 2**32
pi = F(22,7)
error = max(pi**2/(4*n_disk**2), pi**2/(8*k_arc**2),
            F(99,100)*pi**2/(2*n_smooth**2)) + F(12,den)
bound = (minimum - pi*error**2) / (1 + 2*error/(F(1,2)-error))
assert bound > T
sizes = [n_disk] + [n*((k_arc+n-1)//n) for n in (3,5,7,9)] + [n_smooth]
print(json.dumps(dict(target=str(T),continuous_lower=str(minimum),rounding_grid=den,
    denominator=den+10,disk_samples=n_disk,arc_resolution=k_arc,smooth_samples=n_smooth,
    group_size_upper=sizes,total_size_upper=sum(sizes),error=str(error),
    finite_lower=str(bound),positive_margin=str(bound-T)),indent=2))
