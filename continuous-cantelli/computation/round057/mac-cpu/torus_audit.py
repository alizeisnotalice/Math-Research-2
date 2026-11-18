"""Independent saddle-centered Cauchy integrals for asymmetric Wick groups."""
import json,math
from pathlib import Path

import numpy as np
from scipy.special import gammaln,logsumexp
from scipy.optimize import root


def compare(m,N,nu,alpha,beta,x0):
    d=m-1
    k=round(N*nu)
    nu=k/N
    t=2*m*k
    M=N//2-d*k
    lam=d/(2*m-1)*(4*m/(2*m-1))**(m/d)*1.01
    a=N**(-1/(2*d))
    A=abs((2*k*alpha-x0)*a)
    B=(2*k*beta-.5)*a*a
    Rp=math.sqrt((1-2*d*nu)/lam)
    Up=Rp/d
    def stationary(logs):
        R,U=np.exp(logs);s=R+U
        v=t/s+A+2*B*s
        return [(lam*N*R*R+R*v-N)/N,(U*v-2*k)/N]
    ans=root(stationary,np.log([Rp,Up]),tol=1e-11)
    assert max(abs(np.array(stationary(ans.x))))<1e-10
    R,U=np.exp(ans.x);S=R+U
    log_phase=lam*N*R*R/2+t*math.log(S)+A*S+B*S*S-N*math.log(R)-2*k*math.log(U)
    ref=gammaln(t+1)-gammaln(2*k+1)-gammaln(t-2*k+1)+M*math.log(lam*N/2)-gammaln(M+1)
    def cauchy(n):
        theta=np.arange(n)*(2*np.pi/n)
        phase=np.exp(1j*theta)
        W=R*phase[:,None]
        T=U*phase[None,:]
        D=W+T
        with np.errstate(divide='ignore',invalid='ignore'):
            Z=lam*N*(W*W-R*R)/2+t*np.log(D/S)+A*(D-S)+B*(D*D-S*S)
            vals=np.exp(Z-1j*(N*theta[:,None]+2*k*theta[None,:]))
        vals=np.where(D==0,0,vals)
        i=np.mean(vals)
        assert i.real>0 and abs(i.imag/i.real)<1e-8
        return math.log(i.real)+log_phase-ref
    recurrence=[0.,-np.inf]
    weight=0.
    total=0.
    for j in range(1,2*M+1):
        lp=np.logaddexp(math.log(A)+recurrence[0],math.log(2*B)+recurrence[1])-math.log(j)
        recurrence=[lp,recurrence[0]]
        if j%2==0:
            ell=j//2
            weight+=math.log((M-ell+1)*2/(lam*N))+math.log((t+j-1)*(t+j)/((t-2*k+j-1)*(t-2*k+j)))
            total=np.logaddexp(total,lp+weight)
    near=cauchy(512)
    fine=cauchy(1024)
    initial_s=Rp+Up
    C=2*lam*Rp*Rp
    E=2*d*nu/m
    linear=2*nu*abs(alpha)*initial_s
    quadratic=2*nu*beta*initial_s**2-2*nu**2*alpha**2*(initial_s**2/C+Up**2/E)
    return dict(m=m,N=N,k=k,nu=nu,alpha=alpha,beta=beta,x0=x0,A=A,B=B,
                recurrence=float(total),cauchy512=near,cauchy1024=fine,stationary_relative_gain=linear,
                stationary_quadratic_gain=quadratic,scaled_second_order=(total-N*a*linear)/(N*a*a),
                fine_disagreement=float(abs(fine-total)),nested_disagreement=abs(fine-near))


def main():
    rs=[compare(m,N,.43/(2*(m-1)),x,.6,.3) for m in (3,5,8) for N in (128,512,2048,4096)
        for x in (-.55,.37)]
    Path('torus_audit.json').write_text(json.dumps(rs,indent=2)+"\n")
    print(max(x['fine_disagreement'] for x in rs),max(x['nested_disagreement'] for x in rs))

if __name__=='__main__': main()
