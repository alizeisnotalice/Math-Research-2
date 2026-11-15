from fractions import Fraction as Q
import math,sys,json
def U(x):return 1/(abs(x[0])+abs(x[1]))
def sub(x,z):return (x[0]-z[0],x[1]-z[1])
nu=[(Q(1,6),(Q(3),Q(0))),(Q(1,6),(-Q(3),Q(0))),(Q(1,6),(Q(0),Q(3))),(Q(1,6),(Q(0),-Q(3))),(Q(1,6),(Q(2),Q(2))),(Q(1,6),(-Q(2),-Q(2)))]
assert sum(w for w,z in nu)==1
assert all(sum(w*z[j] for w,z in nu)==0 for j in (0,1))
L=Q(3);rows=[]
for x in [(Q(10),Q(12)),(-Q(10),Q(12)),(Q(1,2),Q(15)),(Q(15),-Q(1,2))]:
 rho=abs(x[0])+abs(x[1]);W=U(x)-sum(w*U(sub(x,z)) for w,z in nu)
 general=2*L/(rho*(rho-2*L));assert abs(W)<=general
 corner=abs(x[0])>L and abs(x[1])>L
 if corner:
  bound=4*L*L/(rho*rho*(rho-2*L));assert W<=0 and abs(W)<=bound
 else:bound=general
 rows.append(dict(x=list(map(str,x)),W=str(W),bound=str(bound),corner=corner))
a=Q(3,2);bad=Q(1)-(1/(a-1)+1/(a+1))/2
assert bad==-Q(1,5)
N=8192;quad=[]
for a,b in [(1.,1.),(1.,2.),(.5,3.)]:
 def f(th):
  s=math.sin(th)**2;c=math.cos(th)**2
  return s/((s+a*a*c)*(s+b*b*c))
 h=(math.pi/2)/N;v=(f(0)+f(math.pi/2)+4*sum(f(h*i) for i in range(1,N,2))+2*sum(f(h*i) for i in range(2,N,2)))*h/3
 exact=math.pi/(2*(a+b));assert abs(v-exact)<1e-11
 quad.append(dict(a=a,b=b,quadrature=v,closed_form=exact,abs_error=abs(v-exact)))
print(json.dumps(dict(python=sys.version,scope="exact rational tails and domination counterexample; numeric scalar semigroup integral only",probability_and_zero_mean=True,tail_samples=rows,separated_dipole_W_at_vertex=str(bad),rho_min_3_sufficient_m="1/2",potential_integral_Simpson_N=N,potential_integral_samples=quad),indent=2))
