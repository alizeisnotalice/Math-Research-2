from fractions import Fraction as F
import mpmath as mp
import sympy as S
mp.mp.dps=70
cases={3:['0.73585385','0.79612902','0.97037297'],4:['0.65074430','0.69319210','0.91118522','0.96380609']}
for n,start in cases.items():
 def c(k): return F(int(S.catalan(k)),4**k)
 def mF(*xx):
  j0=(-1)**n+2*sum((-1)**i*xx[i] for i in range(n))
  return tuple(((-1)**n+2*sum((-1)**i*xx[i]**(2*k+1) for i in range(n)))/(2*k+1)-mp.mpf(c(k).numerator)/c(k).denominator*j0 for k in range(1,n+1))
 root=mp.findroot(mF,tuple(map(mp.mpf,start)),tol=mp.mpf('1e-60'),maxsteps=100)
 print('n',n,'root',*[mp.nstr(x,22) for x in root])
 center=[F(str(mp.nstr(x,17))) for x in root]; eps=F(1,10**12)
 j0=(-1)**n+2*sum((-1)**i*center[i] for i in range(n))
 fv=S.Matrix([S.Rational(((-1)**n+2*sum((-1)**i*center[i]**(2*k+1) for i in range(n)))/(2*k+1)-c(k)*j0) for k in range(1,n+1)])
 jac=S.Matrix([[S.Rational(2*(-1)**i*(center[i]**(2*k)-c(k))) for i in range(n)] for k in range(1,n+1)])
 A=jac.inv(); af=A*fv
 # J_ki=2(-1)^i(x_i^(2k)-c_k). Bound deviations from center by derivative max.
 bounds=[[F(4*k)*(center[i]+eps)**(2*k-1)*eps for i in range(n)] for k in range(1,n+1)]
 margin=[]
 for i in range(n):
  correction=F(abs(af[i])); contraction=sum(sum(F(abs(A[i,k]))*bounds[k][j] for k in range(n)) for j in range(n))
  margin.append((correction+eps*contraction)/eps)
 assert max(margin)<1
 assert center[0]>eps and center[-1]<1-eps and all(center[i]+2*eps<center[i+1] for i in range(n-1))
 def atan_bounds(x, terms):
  partial=sum(((-1)**j)*x**(2*j+1)/F(2*j+1) for j in range(terms))
  following=partial+((-1)**terms)*x**(2*terms+1)/F(2*terms+1)
  return min(partial,following),max(partial,following)
 l5,h5=atan_bounds(F(1,5),80); l239,h239=atan_bounds(F(1,239),20)
 pi_lo=16*l5-4*h239; pi_hi=16*h5-4*l239
 lo=pi_lo/(4*(j0+2*n*eps)); hi=pi_hi/(4*(j0-2*n*eps))
 print('center',*[str(x) for x in center], 'eps',eps)
 print('contraction+correction max ratio',float(max(margin)))
 target={3:(F('0.9575742378958'),F('0.9575742379100')),4:(F('0.9697918061584'),F('0.9697918061777'))}[n]
 assert target[0]<lo<hi<target[1]
 print('lambda interval',float(lo),float(hi),'certified target',*[str(x) for x in target])
 print('mF infnorm',mp.nstr(max(abs(v) for v in mF(*root)),5))
