from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import platform

ROOT=Path(__file__).parent


def main():
    pi_min,pi_max=Q(157,50),Q(22,7)
    inv_a_bound=Q(8,3)
    m_min=Q(5,4)
    qmax=Q(1,20000)
    assert 2*pi_max<inv_a_bound**2 and 2*m_min**2<pi_min
    assert Q(8,3)**8>2500
    log_pi_bound=Q(23,20)
    assert sum(log_pi_bound**i/factorial(i) for i in range(7))>pi_max
    z_loss=4*qmax*inv_a_bound + Q(1,2500)
    moment_loss=16*qmax*inv_a_bound+Q(17,4*2500)
    assert z_loss<Q(1,1000) and moment_loss<Q(1,250)
    mass_min=Q(999,1000)
    c_q_min=m_min-Q(1,250)
    d_q_max=log_pi_bound+Q(2,999)
    test=Q(9,5)
    heat_min=mass_min*(test*c_q_min-d_q_max)/(test**2/2)
    free_max=Q(2,3)/(1-qmax)
    margin=heat_min-free_max
    assert margin>Q(59,10000)
    numeric=dict(gaussian_weight_loss_upper=str(z_loss),weighted_first_moment_loss_upper=str(moment_loss),
                 conditional_first_moment_lower=str(c_q_min),entropy_upper_part=str(d_q_max),
                 negative_term_lower_exact=str(heat_min),positive_term_upper_exact=str(free_max),
                 margin_exact=str(margin),margin_display=float(margin),python=platform.python_version())
    (ROOT/'certificate.json').write_text(json.dumps(numeric,indent=2)+'\n')
    print(json.dumps(numeric))


if __name__=='__main__':
    main()
