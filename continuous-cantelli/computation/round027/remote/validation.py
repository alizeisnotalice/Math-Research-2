"""High-precision endpoint and width checks for the tilted plateau envelope."""
import json
from pathlib import Path
import platform

import mpmath as mp


def tail(q):
    return mp.erfc(q/mp.sqrt(2))/2


def constants(a,m,b,r,lam):
    z=mp.sqrt(1+a*lam**2)
    core=mp.exp(lam**2/(2*z*z))/z*(tail(-r*z+lam/z)-tail(r*z+lam/z))
    q_m=mp.exp((1-m)*lam**2/2)
    q_b=mp.exp((1-b)*lam**2/2)
    p_core=tail(lam-r)-tail(lam+r)
    return core+q_b-q_m*p_core, q_m-q_b


def solve_tail_log(log_p, guess):
    g=mp.mpf(guess)
    lo=mp.mpf(0)
    hi=max(mp.mpf(1),2*g)
    assert mp.log(tail(lo))>log_p>mp.log(tail(hi))
    for _ in range(220):
        mid=(lo+hi)/2
        if mp.log(tail(mid))>log_p:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def boundary_integral(z,power):
    # Shifts to outward distance t from the cutoff, avoiding
    # subtraction of large Gaussian exponents.
    b=mp.exp(power-z**2/2)/mp.sqrt(2*mp.pi)
    return b*mp.quad(lambda t:mp.exp(-z*t-t*t/2),[0,1/max(1,z),1,4,mp.inf])


def control_f(a,m,b,r,lam,d1,d2):
    c0,scale=constants(a,m,b,r,lam)
    return c0+scale*(tail(lam*d1)+tail(lam*d2))


def main():
    mp.mp.dps=75
    a,r,b=mp.mpf(1),mp.mpf(1),mp.mpf(2)
    m=mp.mpf(3)/4
    s=mp.sqrt(1-m)
    c_a=mp.exp(1/(2*a))/mp.sqrt(a)
    evaluations=[]
    for lam in [20,32,128,512,2048,8192,32768]:
        lam=mp.mpf(lam)
        c0,scale=constants(a,m,b,r,lam)
        q=solve_tail_log(mp.log1p(-c0)-mp.log(scale)-mp.log(2),s*lam)
        d=q/lam
        approx1=s-mp.log(lam*s*mp.sqrt(2*mp.pi)/2)/(s*lam**2)
        approx2=approx1+c_a/(s*lam**3)
        exponent=(1-m)*lam**2/2
        analytic=mp.exp(exponent)*tail(q)
        integral=boundary_integral(q,exponent)
        evaluations.append(dict(lam=str(lam),d=str(d), p=str(control_f(a,m,b,r,lam,d,d)),
                                scaled_correction=str((d-approx1)*lam**3),
                                scaled_remainder=str((d-approx2)*lam**4/mp.log(lam)**2),
                                scaled_reference=str(c_a/s),
                                independent_relative_error=str(abs(integral/analytic-1))))
    transition=[]
    for beta in [mp.mpf('0.01'),mp.mpf(1),mp.mpf(5),mp.mpf(25)]:
        q0=solve_tail_log(-beta*beta/2-mp.log(2),beta)
        shift=c_a*mp.exp(-beta**2/2+q0*q0/2)*mp.sqrt(2*mp.pi)/2
        for lam in [32,512,8192,32768]:
            lam=mp.mpf(lam)
            c0,scale=constants(a,1-(beta/lam)**2,b,r,lam)
            q=solve_tail_log(mp.log1p(-c0)-mp.log(scale)-mp.log(2),q0)
            transition.append(dict(beta=str(beta),lam=str(lam),q=str(q),q0=str(q0),
                                   scaled_shift=str((q-q0)*lam), reference_shift=str(shift)))
    lam=mp.mpf(20)
    cuts={str(d):str(control_f(a,m,b,r,lam,d,d)) for d in (mp.mpf('0.485'),mp.mpf('0.49'))}
    out=dict(python=platform.python_version(), mpmath=mp.__version__, precision_digits=mp.mp.dps,
             fixed_floor=evaluations, near_one_transition=transition, finite_cuts=cuts)
    Path("validation.json").write_text(json.dumps(out,indent=2)+"\n")
    print(cuts, evaluations[0]["d"], evaluations[-1]["scaled_correction"], str(c_a/s))


if __name__=='__main__':
    main()
