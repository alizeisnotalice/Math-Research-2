"""Exact illustrative finite-cell C4 certificate; no optimized upper bound."""
import hashlib,json,pathlib,sys
from fractions import Fraction as Q

def atan_bounds(q,n):
    terms=[(-1)**k*q**(2*k+1)/Q(2*k+1) for k in range(n)]
    s=sum(terms,Q());nxt=(-1)**n*q**(2*n+1)/Q(2*n+1)
    return min(s,s+nxt),max(s,s+nxt)
a,b=atan_bounds(Q(1,5),10);c,d=atan_bounds(Q(1,239),4)
alpha,beta=16*a-4*d,16*b-4*c
assert Q(3)<alpha<beta<Q(22,7)
h=Q(1);B=Q(2);N=4
cells=[((Q(a),Q(b)),(Q(c),Q(d))) for a,b in [(-2,0),(0,2)] for c,d in [(-2,0),(0,2)]]
areas=[(r[0][1]-r[0][0])*(r[1][1]-r[1][0]) for r in cells]
def maxabs(interval):return max(map(abs,interval))
def diffmax(I,J):return max(abs(x-y) for x in I for y in J)
def lower_density(d1,d2):return h*h/(beta*beta*(h*h+d1*d1)*(h*h+d2*d2))
bvec=[lower_density(maxabs(r[0]),maxabs(r[1])) for r in cells]
L=[[areas[j]*lower_density(diffmax(r[0],s[0]),diffmax(r[1],s[1])) for j,s in enumerate(cells)] for r in cells]
def matvec(v):return [sum((q*w for q,w in zip(row,v)),Q()) for row in L]
a=[Q() for _ in cells];term=bvec
for _ in range(N):
    a=[x+y for x,y in zip(a,term)]
    term=matvec(term)
La=matvec(a);residual=[bvec[i]+La[i]-a[i] for i in range(len(a))]
assert residual==term and all(x>0 for x in residual)
assert all(a[i]==a[len(a)-1-i] for i in range(len(a)))
q=2/beta*(h/(2*B))/(1+(h/(2*B))**2)
assert max(map(sum,L))<=1-q
mu_mass=1+sum((v*area for v,area in zip(a,areas)),Q())
Qstr=lambda x:str(x)
out=dict(python=sys.version,scope='four-cell exact source-positivity/C4 illustration; not optimized and no D minimum certification',script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),h=str(h),Omega='[-2,2]^2',Neumann_terms=N,pi_lower=str(alpha),pi_upper=str(beta),pi_width_float=float(beta-alpha),cells=[[[str(x) for x in interval] for interval in r] for r in cells],b=list(map(Qstr,bvec)),L=[list(map(Qstr,row)) for row in L],a=list(map(Qstr,a)),residual_b_plus_La_minus_a=list(map(Qstr,residual)),a_float=list(map(float,a)),residual_float=list(map(float,residual)),row_sum_bound_float=float(max(map(sum,L))),survival_q_lower_float=float(q),mu_mass_exact=str(mu_mass),mu_mass_float=float(mu_mass),kernel_integral_formula='2*pi*h*mu_mass',kernel_half_integral_over_pi_exact=str(h*mu_mass),source_formula='nu=p_h+p_h*f-f >=0; integral nu=1; f=sum a_i 1_Ci')
print(json.dumps(out,indent=2))
