from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json


def exp_bound(x,n=36):
    terms=[x**j/factorial(j) for j in range(n+2)]
    lower=sum(terms[:-1])
    upper=lower+terms[-1]/(1-x/(n+2))
    return lower,upper


def gaussian_integral(n):
    return sum(Q((-1)**j,2**j*factorial(j)*(2*j+1)) for j in range(n+1))


def main():
    pi_lo,pi_hi=Q(157,50),Q(22,7)
    sigma_lo,sigma_hi=Q(5,2),Q(251,100)
    assert sigma_lo*sigma_lo<2*pi_lo<2*pi_hi<sigma_hi*sigma_hi
    e0_lo,e0_hi=exp_bound(Q(1,2))
    phi_lo,phi_hi=1/(sigma_hi*e0_hi),1/(sigma_lo*e0_lo)
    m_lo=1-2*gaussian_integral(12)/sigma_lo
    m_hi=1-2*gaussian_integral(13)/sigma_hi

    t_lo,t_hi=pi_lo*pi_lo/8,pi_hi*pi_hi/8
    e_lo=1/exp_bound(t_hi)[1]
    e_hi=1/exp_bound(t_lo)[0]
    # Odd modes k=3,5,... have pairwise square differences at least 16.
    # Their excitation ratio is at most 1/2 by exp(x) >= 1 + x.
    assert 1+16*t_lo>2
    tail=2/exp_bound(9*t_lo)[0]
    capmass_lo=4*e_lo/pi_hi-2*tail
    j_hi=(4/pi_lo-32/pi_hi**3)*e_hi+Q(2,3)*tail
    previous_upper=2*phi_hi-m_lo
    required_lower=2*phi_lo-j_hi
    available_lower=1-2*m_hi
    available_upper=1-2*m_lo

    assert capmass_lo>m_hi and previous_upper<available_lower
    assert required_lower>available_upper

    def lower(x):
        return x.numerator*1000000//x.denominator

    intervals={
        'cap_mass_lower':capmass_lo,
        'required_mass_upper':m_hi,
        'R84_budget_upper':previous_upper,
        'budget_available_lower':available_lower,
        'budget_available_upper':available_upper,
        'R88_budget_lower':required_lower,
        'positive_exclusion_margin':required_lower-available_upper,
    }
    decimal_bounds={}
    for label,num in intervals.items():
        outward=lower(num)+(label.endswith('_upper'))
        decimal_bounds[label]=f'{outward//1000000}.{outward%1000000:06d}'
    report={'A':1,'B':1,'T':1,'outward_decimal_bounds':decimal_bounds,
            'all_checks_passed':True}
    Path(__file__).with_name('exact_example.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':
    main()
