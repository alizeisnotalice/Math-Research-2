"""Calibrate finite well identities and illustrative critical well budgets."""

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp


def exact_algebra():
    t,s,x,f = sp.symbols('t s x f')
    odds = (t-s)*x-(t*t-s*s)*f/2
    expected = (t-s)*(x-t*f)+(t-s)**2*f/2
    assert sp.expand(odds-expected)==0
    logp = t*x-t*t*f/2
    assert sp.diff(logp,t)==x-t*f and sp.diff(logp,t,2)==-f
    return str(sp.expand(odds))


def gauss_partition(dps, method):
    with mp.workdps(dps):
        cut = (-mp.inf, mp.mpf('-1'), mp.mpf('.3'), mp.mpf('2'), mp.inf)
        t,u = mp.mpf('1.3'),mp.mpf('-.35')
        g = lambda x: mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
        Phi = lambda x: mp.erfc(-x/mp.sqrt(2))/2
        rows = []
        cats = mp.mpf(0)
        sum_score = mp.mpf(0)
        sum_variance = mp.mpf(0)
        maxdiff = mp.mpf(0)
        for lo,hi in zip(cut[:-1],cut[1:]):
            p=Phi(hi-t)-Phi(lo-t)
            q=Phi(hi-u)-Phi(lo-u)
            score=(g(lo-t)-g(hi-t))/p
            def yphi(y):
                return y*g(y) if mp.isfinite(y) else mp.mpf(0)
            m2=1+(yphi(lo-t)-yphi(hi-t))/p
            variance=m2-score**2
            formula=(t-u)*score+(t-u)**2/2+mp.log(q/p)
            direct=mp.quad(lambda x: g(x-t)/p*((t-u)*x-(t*t-u*u)/2+mp.log(q/p)),
                           (lo,hi),method=method)
            maxdiff=max(maxdiff,abs(formula-direct))
            cats+=p*mp.log(p/q)
            sum_score+=p*score
            sum_variance+=p*m2
            rows.append((p,q,score,variance,formula,direct))
        assert abs(sum_score)<mp.mpf(10)**(-dps+4)
        assert abs(sum_variance-1)<mp.mpf(10)**(-dps+4)
        assert maxdiff<mp.mpf(10)**(-dps+6)
        cond=sum(row[0]*row[4] for row in rows)
        assert abs(cats+cond-(t-u)**2/2)<mp.mpf(10)**(-dps+5)
        return rows,cats,cond,maxdiff


def illustrative():
    L=mp.mpf('10000')
    A=2*L
    K=mp.sqrt(mp.log(A))
    zq=mp.sqrt(2*mp.pi)*mp.erf(K/mp.sqrt(2))+2*mp.exp(-K*K/2)/K
    Phi = lambda x: mp.erfc(-x/mp.sqrt(2))/2
    vals=[]
    for a in (mp.mpf('.5'),mp.mpf('.8'),mp.mpf('.95')):
        sd=mp.sqrt(1-a*a)
        windows=((-a-sd,-a+sd),(a-sd,a+sd))
        if windows[0][1]>windows[1][0]:
            windows=((windows[0][0],windows[1][1]),)
        p=mp.mpf(0)
        q=mp.mpf(0)
        for lo,hi in windows:
            assert max(abs(lo),abs(hi))<K
            p+=(Phi((hi-a)/sd)-Phi((lo-a)/sd)
                +Phi((hi+a)/sd)-Phi((lo+a)/sd))/2
            q+=mp.sqrt(2*mp.pi)/zq*(Phi(hi)-Phi(lo))
        class_kl=p*mp.log(p/q)+(1-p)*mp.log((1-p)/(1-q))
        delta=1-mp.exp(-a*a/2)*mp.cosh(a)
        vals.append({
            'symmetric_mixture_base_a':mp.nstr(a,10), 'union_probability':mp.nstr(p,35),
            'Huber_probability':mp.nstr(q,35), 'coarse_entropy_required_delta':mp.nstr(2*class_kl,35),
            'local_tilt_defect_at_lambda_one':mp.nstr(delta,35),
            'indicative_required_energy_ratio_for_log_lr_threshold_100':
                mp.nstr(delta/(mp.mpf('4e-4')+mp.mpf('1e-2')),35),
        })
    return vals


def main():
    algebra = exact_algebra()
    primary=gauss_partition(66,'tanh-sinh')
    alternate=gauss_partition(48,'gauss-legendre')
    with mp.workdps(66):
        cross=max(abs(a-b) for row,other in zip(primary[0],alternate[0]) for a,b in zip(row,other))
        assert cross<mp.mpf('1e-35'),cross
        keys=('p_t','p_u','score_mean','score_variance','conditional_kl_formula','conditional_kl_integral')
        output={'scope':'Exact Gaussian tilt benchmark and illustrative symmetric mixture profiles; the mixtures are not exact graphs.',
                'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,
                'likelihood_algebra':algebra,'true_gaussian_cells':
                    [{k:mp.nstr(v,37) for k,v in zip(keys,row)} for row in primary[0]],
                'categorical_kl':mp.nstr(primary[1],38),'mean_conditional_kl':mp.nstr(primary[2],38),
                'max_cross_method_difference':mp.nstr(cross,14),'illustrative_multiwell':illustrative()}
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('multiwell-results.json')
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(output['max_cross_method_difference'],output['categorical_kl'],output['mean_conditional_kl'])


if __name__=='__main__':
    main()
