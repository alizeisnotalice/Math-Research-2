"""Combine the independent AH and EH area gains by exact arithmetic."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path


def run(ah_path, eh_path):
    ah, eh = [json.loads(p.read_bytes()) for p in (ah_path,eh_path)]
    assert ah['complete'] and eh['complete']
    assert ah['baseline_sha256']==eh['baseline_round222_sha256']
    delta0,delta1 = map(F,ah['common_AH_angle_band'])
    assert F(eh['EH_AH_separation']['exact'][0]) > F('.45549917') > F(6,5)*(delta1-delta0)
    gain_lo, gain_hi = map(F,ah['deeper_AH_exact_gain']['exact'])
    assert 0 < gain_lo < gain_hi
    eh_lo, eh_hi = map(F,eh['area']['most_cut_corner']['exact'])
    lower, upper = eh_lo-gain_hi, eh_hi-gain_lo

    def outward(a, round_up):
        k = 10**42
        d = (a*k).__ceil__() if round_up else (a*k).__floor__()
        return f'{d//k}.{d%k:042d}'

    return dict(complete=True, AH_band=ah['common_AH_angle_band'], EH_low_corner=eh['low_corner'],
                EH_high_corner=eh['high_corner'],
                product_wide_area_upper=eh['area']['least_cut_corner']['decimal'][1],
                most_cut_corner_area=dict(exact=list(map(str,(lower,upper))),
                                          decimal=[outward(lower,False),outward(upper,True)]),
                input_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in (ah_path,eh_path)})


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    for name in ('ah','eh','output'):
        p.add_argument(name,type=Path)
    a = p.parse_args()
    result = run(a.ah,a.eh)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['most_cut_corner_area']['decimal'])
