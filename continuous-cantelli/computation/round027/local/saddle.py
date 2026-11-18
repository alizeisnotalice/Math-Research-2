"""Independent high-precision real and complex critical saddle checks."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json, math, os
from pathlib import Path
import platform, re, subprocess, sys, time

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import qmc

R = .25
SEED = 251118
S_EXPONENTS = (4, 6, 8, 10, 12, 14)
WORKERS = 12
CPU_RE = re.compile(r"CPU usage: (.*?)% user, (.*?)% sys, (.*?)% idle")
DPS = 58


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def profiles(size):
    grid = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2((size-1).bit_length())[:size]
    for i, (hv, bv, cv, xv) in enumerate(grid):
        n = (4, 6, 8)[i % 3]
        h = 0 if i % 4 == 0 else float(2.5*hv-1.25)
        b = float(3*bv-1.5)
        c = float(math.exp(math.log(4)*cv-math.log(2)))
        x = float(6*xv-2)
        yield i, n, h, b, c, x


def calc_constants(n):
    k = n-2
    theta = mp.pi/(2*(n-1))
    analytic = (n-1)*mp.sin(theta)/n
    root = (k*mp.sin(theta)/2)**(mp.mpf(1)/n)
    D = root**(n+k)/k
    check = mp.mpf(2*(n-1))/(n*k)*(k*D)**(mp.mpf(n)/(2*n-2))
    assert mp.almosteq(check, analytic)
    return k, theta, analytic, root, D


def real_mp(n, k, D, root, lam, s, H, b, x):
    coeff = (D, mp.mpf(1)/n, H*s/n, b*s*s/n)
    powers = (-k, n, n+1, n+2)

    def derivative(z, j):
        return mp.fsum(c * math.prod(p-i for i in range(j)) * z**(p-j)
                       for c,p in zip(coeff,powers))

    m = mp.findroot(lambda z: derivative(z, 1), (root*.95, root*1.05))
    assert 0 < m < R/s
    curvature = derivative(m, 2)
    assert curvature > 0
    width = 1/mp.sqrt(lam*curvature)
    minimum = derivative(m, 0)
    lo = max(-22, -m/width)
    hi = min(36, (R/s-m)/width)
    cuts = [lo] + [j for j in (-11, -4, 0, 4, 11, 22) if lo < j < hi] + [hi]

    def sq(y):
        return lam*(derivative(m+width*y,0)-minimum)

    def integrand(y):
        u = s*(m+width*y)
        return mp.exp(-sq(y)-x*u-u*u/2) if u > 0 else mp.mpf(0)

    val = mp.quad(integrand, cuts)
    tail = mp.mpf(0)
    for y, direction, clipped in [(lo,-1,lo == -22),(hi,1,hi == 36)]:
        if clipped:
            ay = mp.exp(x*x/2-sq(y))
            dy = direction*lam*width*derivative(m+width*y,1)
            assert dy > 0
            tail += ay/dy
    return m, curvature, minimum, width, val, tail


def fast_real(n,k,D,m,lam,s,H,b,x,width):
    """Evaluate S(y)=lambda*(E(m+width*y)-E(m)) with exact linear cancellation."""
    power_terms = ((n, m**n/n), (n+1,H*s*m**(n+1)/n), (n+2,b*s*s*m**(n+2)/n))
    neg_poly = tuple(k*math.comb(k,j-1)-(math.comb(k,j) if j<=k else 0)
                     for j in range(2,k+2))
    polys = [(coef, tuple(math.comb(p,j) for j in range(2,p+1))) for p,coef in power_terms]

    def pol(coeffs, z):
        total = 0
        for coef in reversed(coeffs):
            total = total*z+coef
        return total

    lo = max(-22, -m/width)
    hi = min(36, (R/s-m)/width)

    def integrand(y):
        z = width*y/m
        if z <= -1:
            return 0.
        remainder = D*m**(-k)*pol(neg_poly,z)/(1+z)**k
        remainder += math.fsum(coef*pol(poly,z) for coef,poly in polys)
        action = lam*z*z*remainder
        u = s*(m+width*y)
        return math.exp(-action-x*u-u*u/2)

    points = [t for t in (-11,-4,0,4,11,22) if lo<t<hi]
    return quad(integrand,lo,hi,points=points,epsabs=1e-12,epsrel=1e-12,limit=80)[0]


def complex_mp(n,theta,s,h,b):
    values = []
    for initial in (mp.exp(1j*theta), -mp.exp(-1j*theta)):
        def diff(z):
            return 1j - z**(n-1) - (n+1)*h*s*z**n/n - (n+2)*b*s*s*z**(n+1)/n
        z = mp.findroot(diff, (initial, initial*(1+mp.mpf(".01"))))
        action = 1j*z - (z**n)*(1+h*s*z+b*s*s*z*z)/n
        second = -(n-1)*z**(n-2) - (n+1)*h*s*z**(n-1) - (n+2)*(n+1)*b*s*s*z**n/n
        values.append((z, action, second))
    return values


def case_rows(case, smoke):
    index, n, h, b, c, x = case
    h, b, c, x = map(mp.mpf, (h,b,c,x))
    k, theta, A, root, D = calc_constants(n)
    # v~c*u^n, xi=2/(n*c*s^(n-1)), lambda=xi*s.
    linear = abs(h)*(mp.sin(2*theta)-root**(n+1))/n
    A2 = (b/n-(n+1)**2*h*h/(2*n*n*(n-1)))*mp.sin(3*theta)
    E2 = (b/n-(n+1)**2*h*h/(4*n*n*(n-1)))*root**(n+2)
    quadratic = A2+E2
    rows = []
    for ep in ((4,14) if smoke else S_EXPONENTS):
        s = mp.mpf(2)**(-ep)
        lam = 2/(n*c*s**k)
        halves = [real_mp(n,k,D,root,lam,s,sigma*h,b,sigma*x) for sigma in (1,-1)]
        complex_sides = complex_mp(n,theta,s,h,b)
        E = min(side[2] for side in halves)
        F = max(mp.re(side[1]) for side in complex_sides)
        reals = []
        error_max = 0
        for sigma,(m,curv,minimum,width,integral,tail) in zip((1,-1),halves):
            faster = fast_real(n,k,*map(float,(D,m,lam,s,sigma*h,b,sigma*x,width)))
            error = abs(float(integral)-faster)
            error_max = max(error_max,error)
            if error > 2e-10*max(1,abs(faster)) or tail > mp.mpf("2e-10"):
                raise ArithmeticError((case,ep,sigma,error,str(tail)))
            ratio = integral / (mp.sqrt(2*mp.pi)*mp.exp(-sigma*x*s*m-(s*m)**2/2))
            reals.append({
                "sign":sigma,"m":mp.nstr(m,24),"min":mp.nstr(minimum,24),
                "curvature":mp.nstr(curv,24),"normalized_integral":mp.nstr(integral,24),
                "integral_to_local_gaussian":mp.nstr(ratio,24),
                "tail_after_normalizing_bound":mp.nstr(tail,12),
            })
        rows.append({
            "s_exponent":-ep,"lambda":float(lam),"real":reals,
            "complex_real": [mp.nstr(mp.re(v[1]),24) for v in complex_sides],
            "complex_curvature_real": [mp.nstr(mp.re(-v[2]),24) for v in complex_sides],
            "advantage":mp.nstr(F+E,24),"advantage_times_lambda":mp.nstr(lam*(F+E),24),
            "linear_error_scaled":float(abs((F+E)/s-linear)),
            "quadratic_error_scaled":float(abs((F+E-linear*s)/s**2-quadratic)),
            "cross_error":error_max,
        })
    return {"index":index,"n":n,"h":float(h),"b":float(b),"c":float(c),"x":float(x),
            "D":mp.nstr(D,24),"d":mp.nstr(2*D/(n*c),24),
            "leading_action":mp.nstr(A,24),"linear_gap":mp.nstr(linear,24),
            "quadratic_gap":mp.nstr(quadratic,24),"rows":rows}


def run_shard(batch,smoke):
    mp.mp.dps = DPS
    start = utc()
    cpu = time.process_time()
    rows = [case_rows(case,smoke) for case in batch]
    return {"pid":os.getpid(),"start":start,"finish":utc(),"cpu":time.process_time()-cpu,"rows":rows}


def snapshot(pid):
    process = subprocess.check_output(["/bin/ps","-axo","pid=,ppid=,%cpu=,time=,state="],text=True)
    return [row for line in process.splitlines() if (row:=line.split())[1]==str(pid)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases",type=int,default=512)
    parser.add_argument("--out",default="saddle-results.json")
    parser.add_argument("--smoke",action="store_true")
    args = parser.parse_args()
    grid = list(profiles(args.cases))
    if args.smoke:
        grid = grid[:4]
    workers = 4 if args.smoke else WORKERS
    pid = os.getpid()
    started = utc()
    wall = time.monotonic()
    print(json.dumps({"pid":pid,"start":started,"workers":workers,"cases":len(grid)}),flush=True)
    samples = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        fs = [pool.submit(run_shard,grid[j::workers],args.smoke) for j in range(workers)]
        while not all(f.done() for f in fs):
            top = subprocess.check_output(["/usr/bin/top","-l","2","-s","1","-n","0"],text=True)
            sample = {"utc":utc(),"elapsed":time.monotonic()-wall,
                      "cpu_pct_user_sys_idle":list(map(float,CPU_RE.findall(top)[-1])),
                      "ps_pid_ppid_cpu_time_state":snapshot(pid)}
            samples.append(sample)
            print(json.dumps(sample),flush=True)
            time.sleep(.1)
        results = [f.result() for f in fs]
    finish = utc()
    workpids = {f["pid"] for f in results}
    still_running = set(map(int,subprocess.check_output(["/bin/ps","-axo","pid="],text=True).split()))
    assert len(workpids)==workers and not workpids.intersection(still_running)
    out = {"start":started,"finish":finish,"wall_s":time.monotonic()-wall,
           "shards":results,"samples":samples,"workers_gone":True,
           "radius":R,"seed":SEED,"mpmath_dps":DPS,"python":sys.version,
           "platform":platform.platform(),"mpmath":mp.__version__,"scipy":scipy.__version__,"numpy":np.__version__}
    Path(args.out).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"finish":finish,"wall_s":out["wall_s"],"workers_gone":True}),flush=True)


if __name__=="__main__":
    main()
