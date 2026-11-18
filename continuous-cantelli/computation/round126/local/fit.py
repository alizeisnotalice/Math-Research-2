#!/usr/bin/env python3
"""Constrained finite compensation in nested piecewise-linear graph classes."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.optimize import Bounds, LinearConstraint, minimize

NORM = 1 / np.sqrt(2 * np.pi)
CASES = [(1.,.25,.25),(1.,.25,.125),(.75,.25,.25),(1.5,.25,.25)]
LAYERS = [
    (2, (.5, 1.), .18),
    (3, (.75, 1.5, 2.), .12),
    (4, (.375, 1.125, 1.75, 2.5, 3.), .08),
    (6, (.625, .875, 1.375, 2.25, 3.5, 4.), .055),
    (8, (1.625, 2.75, 3.25, 4.5, 5., 6.), .045),
    (10, (.25, 1.875, 3.75, 5.5, 7., 8.), .035),
]


def quadrature(grid, order):
    q, w = leggauss(order)
    z, w = (q + 1) / 2, w / 2
    dx = grid[1] - grid[0]
    x = (grid[:-1,None] + dx*z).ravel()
    basis = np.zeros((len(x),len(grid)))
    i = np.arange(len(x))
    basis[i, i//order] = np.tile(1-z, len(grid)-1)
    basis[i, i//order+1] = np.tile(z, len(grid)-1)
    weight = np.tile(w*dx*NORM, len(grid)-1)
    return x, weight, basis


def transform(parameters, x, weight, basis, b):
    z = np.asarray(parameters, complex)
    p = z[:,None].real
    q = z[:,None].imag
    diff = x[None,:] - p
    # Return the residual multiplied by exp(-Im(z)^2/2).
    kernel = weight[None,:] * np.exp(-diff**2/2 + 1j*q*diff)
    u = -z[:,None]**2/2
    return np.sum(kernel*np.expm1(u*(basis@b)), axis=1)


def solve(grid, previous, params, lip, step):
    x, mass, basis = quadrature(grid, 12)
    t = np.array(params)
    u = -t[:,None]**2 / 2
    weight = mass[None,:] * np.exp(-(x[None,:]-t[:,None])**2/2) / np.minimum(-u,1.)
    lower, upper = np.maximum(-.5,previous-step), np.minimum(.5,previous+step)
    for anchor, value in ((grid[0],0),(grid[-1],0),(0,-.25),(1,0)):
        ix = np.flatnonzero(grid == anchor)[0]
        lower[ix] = upper[ix] = value
    d = np.diff(np.eye(len(grid)), axis=0)
    linear = LinearConstraint(np.c_[d,np.zeros(len(d))], -lip*(grid[1]-grid[0]),
                             lip*(grid[1]-grid[0]))

    def both(v):
        bw = basis @ v[:-1]
        m = np.expm1(u*bw)
        f = np.sum(weight*m,axis=1)
        g = (weight*u*(1+m)) @ basis
        return f,g

    def inequality(v):
        f,_ = both(v)
        return np.r_[v[-1]-f,v[-1]+f]

    def jac(v):
        _,g = both(v)
        return np.c_[np.r_[-g,g],np.ones(2*len(g))]

    eps0 = max(abs(both(np.r_[previous,0.0])[0])) + 1e-8
    objective = lambda v:v[-1]+1e-8*np.sum((v[:-1]-previous)**2)
    grad = lambda v:np.r_[2e-8*(v[:-1]-previous), 1.]
    sol = minimize(objective,np.r_[previous,eps0],jac=grad,method="SLSQP",
                   constraints=(linear,{"type":"ineq","fun":inequality,"jac":jac}),
                   bounds=Bounds(np.r_[lower,0],np.r_[upper,2]),
                   options=dict(maxiter=170,ftol=5e-12))
    v=sol.x
    maxslope=np.max(np.abs(np.diff(v[:-1])))/(grid[1]-grid[0])
    assert maxslope<=lip+2e-10 and max(abs(v[:-1]-previous))<=step+2e-10
    return v,sol,maxslope


def main():
    rows=[]
    for case,(lip,delta,dx) in enumerate(CASES):
        prev_grid=np.arange(-2.,2.+dx/2,dx)
        prev=-delta*np.maximum(0,1-np.abs(prev_grid)/.5)
        ts=[]
        for k,(R,new,step) in enumerate(LAYERS[:int(sys.argv[1])]):
            grid=np.arange(-R,R+dx/2,dx)
            previous=np.interp(grid,prev_grid,prev,left=0,right=0)
            old=ts[:]
            ts=old+list(new)+[-s for s in new]
            v,sol,slope=solve(grid,previous,ts,lip,step)
            b=v[:-1]
            x,w,basis=quadrature(grid,32)
            residual=transform(ts,x,w,basis,b).real
            previous_residual=transform(ts,x,w,basis,previous).real
            row=dict(case=case,layer=k+1,R=R,dx=dx,common_lip=lip,pinned=[0,1,-delta,0],
                     training=ts, training_resid=residual.tolist(),minimax=float(v[-1]),
                     max_abs_training=float(max(abs(residual))),max_abs_previous_change=float(max(abs((residual-previous_residual)[:len(old)]),default=0)),
                     slope=float(slope), l_inf_step=float(max(abs(b-previous))),step_budget=step,
                     scipy_success=bool(sol.success),message=sol.message,iterations=sol.nit,
                     offset=b.tolist())
            rows.append(row)
            print(case,k+1,R,len(ts),row["max_abs_training"],row["max_abs_previous_change"],sol.success,flush=True)
            prev_grid,prev=grid,b
    report=dict(candidates=rows, python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                notes="uniform grid, exterior f=1, derivative cap by adjacent knot differences")
    Path("candidates.json").write_text(json.dumps(report,indent=2)+"\n")
    sizes = [len(row["offset"]) for row in rows]
    max_knots=max(sizes)
    with Path("profiles.bin").open("wb") as f:
        np.array([len(rows), max_knots],dtype="<u4").tofile(f)
        for row in rows:
            values=[row["R"],row["dx"],len(row["offset"]),row["case"],row["layer"]]
            values += row["offset"]+[0.]*(max_knots-len(row["offset"]))
            np.asarray(values,dtype="<f4").tofile(f)


if __name__=="__main__":
    main()
