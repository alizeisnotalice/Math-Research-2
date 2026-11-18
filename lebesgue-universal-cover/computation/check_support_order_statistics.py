"""Numerical cross-checks of the symbolic pair, median and translation formulas."""
from itertools import combinations
import json
import numpy as np
from scipy.spatial import ConvexHull

N = 2**17
theta = np.arange(N) * (2*np.pi/N)
u = np.column_stack((np.cos(theta), np.sin(theta)))
v = np.column_stack((-np.sin(theta), np.cos(theta)))
area = lambda p: ConvexHull(p).volume
h, deriv, boundary = [], [], []
for n, amplitude, phase, shift in [(3,.045,.13,(.02,.01)),(5,.012,.32,(-.01,.03)),(7,.006,.74,(.03,-.02))]:
    assert (n*n-1)*amplitude < .5
    p = .5+amplitude*np.cos(n*(theta-phase))+u@shift
    dp = -n*amplitude*np.sin(n*(theta-phase))+v@shift
    h.append(p)
    deriv.append(dp)
    boundary.append(p[:,None]*u+dp[:,None]*v)

singles = list(map(area,boundary))
pairs, errors = [], []
for i,j in combinations(range(3),2):
    exact_form = (singles[i]+singles[j])/2+np.pi/2*np.mean(abs(h[i]-h[j]))
    pairs.append(area(np.vstack((boundary[i],boundary[j]))))
    errors.append(abs(pairs[-1]-exact_form))
idx = np.argsort(h,axis=0)[1]
m = np.array(h)[idx,np.arange(N)]-.5
dm = np.array(deriv)[idx,np.arange(N)]
rest = np.pi/2*np.mean(dm**2-m**2)
observed = area(np.vstack(boundary))-sum(pairs)/2+np.pi/8
assert max(errors)<1e-7 and abs(rest-observed)<1e-6 and rest>0

checks = 0
for n in (3,5,7,9):
    a = np.pi/(2*n)
    R = .5/np.cos(a)
    # Each sector of angular width 2a is alternately a vertex and its antipodal arc.
    which = np.floor((theta+a)/(2*a)).astype(int)
    centered = (-1.)**which*(R*np.cos(theta-2*a*which)-.5)
    for d in (.001,.03,.07,.2):
        for phase in (0.,.1,.4):
            gain = .5*np.pi*np.mean(abs(centered+d*np.cos(theta-phase))-abs(centered))
            floor = n*d*d/(4*max(.5*np.tan(a)+d,d/a))
            assert gain+1e-9 >= floor
            checks += 1
print(json.dumps(dict(samples=N, maximum_pair_error=max(errors),median_remainder=observed,
                      median_quadrature_error=abs(observed-rest),radial_cases=checks),indent=2))
