"""Exact three-layer relaxation and single-tilt Gaussian calibration (not exact graphs)."""
import json
import platform
import mpmath as mp

mp.mp.dps=80
eps=mp.mpf(3)/5*mp.exp(-mp.mpf(16)/3)

def laplace(t,B):
    # joint expectation exp(t X - B f/2) in the verified three-layer relaxation.
    pieces=[(1,0,1,1),(eps,2,mp.mpf(1)/8,mp.mpf(1)/2),(-eps,2,mp.mpf(1)/8,1),(-eps,2,mp.mpf(5)/8,1),(eps,2,mp.mpf(5)/8,mp.mpf(3)/2)]
    return mp.fsum(weight*mp.exp(mu*t+var*t*t/2-B*time/2) for weight,mu,var,time in pieces)

rows=[]
for t,s,alpha in [(1,2,'0.4'),(2,1,2),(0,1,2),(1,mp.sqrt(2),2),(-2,0,3)]:
    t,s,alpha=map(mp.mpf,(t,s,alpha)); u=alpha*t+(1-alpha)*s; B=alpha*t*t+(1-alpha)*s*s
    gaussian=mp.exp(alpha*(alpha-1)*(t-s)**2/2)
    ratio=laplace(u,B)/gaussian
    rows.append(dict(t=mp.nstr(t,11),s=mp.nstr(s,11),alpha=mp.nstr(alpha,5),B=mp.nstr(B,11),raw_ratio_minus_one=mp.nstr(ratio-1,15),predicted_sign=int(mp.sign(B*alpha*(1-alpha))) if abs(B)>mp.mpf('1e-70') else 0))

floor=[]
for L in (10,30,100,1000):
    v=mp.exp(-2*L)
    def eq(A,B):
        return A-2*L+mp.log(B)+1-(A*A+B)*v,A-B*(1-2*A*v)**2-2*(1-B*v)**2
    A,B=mp.findroot(eq,(2*L-mp.log(2*L)-1,2*L-mp.log(2*L)-3))
    Tmin=A-(1-B*v)-B*(1-A*v)**2/(1-B*v)
    gap=A+mp.log(A)-2*L+1
    floor.append(dict(L=L,A=mp.nstr(A,20),B=mp.nstr(B,20),minimum_T=mp.nstr(Tmin,12),gap=mp.nstr(gap,15),max_moment_residual=mp.nstr(max(map(abs,eq(A,B))),6)))

print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,decimal_precision=mp.mp.dps,renyi_relaxation=rows,gaussian_single_tilt=floor),indent=2))
