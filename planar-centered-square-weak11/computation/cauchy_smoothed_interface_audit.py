"""Scalar checks only; analytic C4/L1 and ray claims are in the proof addendum."""
import hashlib, json, math, pathlib, sys
pi=math.pi

def V(a,b,t):
    a,b=abs(a),abs(b)
    if a==b:
        if a==0: return 2/(pi*t)
        return (math.atan(a/t)+a*t/(t*t+a*a))/(pi*a)
    return 2*(b*math.atan(b/t)-a*math.atan(a/t))/(pi*(b*b-a*a))

def integral(a,b,t,n=8192):
    def f(u):
        q=1-u
        return t**3/((t*t+a*a*q*q)*(t*t+b*b*q*q))
    total=f(0)+f(1)
    total+=4*sum(f(k/n) for k in range(1,n,2))
    total+=2*sum(f(k/n) for k in range(2,n,2))
    return 2*total/(3*n*pi)

rows=[]
for a,b,t in [(0,0,.5),(0,2,.5),(2,2,.5),(.5,2,.75),(3,1,2),(1,1,10)]:
    closed=V(a,b,t);quad=integral(a,b,t)
    assert abs(closed-quad)<1e-11
    rows.append(dict(a=a,b=b,t=t,closed=closed,independent_Simpson=quad,error=abs(closed-quad)))

derivative_rows=[]
for a,b,t in [(.1,.5,.1),(.5,.5,.2),(3,.1,1),(1,0,.5),(1,2,10)]:
    h=1e-5
    d=(V(a+h,b,t)-V(a-h,b,t))/(2*h)
    lower=-1/(a+b)**2
    assert lower-1e-7<=d<=1e-7
    derivative_rows.append(dict(a=a,b=b,t=t,dV_da=d,analytic_lower=lower))

counter_t=1/(10*pi)
interior=2-.5*(V(0,0,counter_t)+V(1,0,counter_t))
assert interior < -8
large=[]
for t in [10.,20.,40.]:
    a,b=1.,2.
    approximation=2/(pi*t)-2*(a*a+b*b)/(3*pi*t**3)
    error=V(a,b,t)-approximation
    large.append(dict(t=t,value=V(a,b,t),cubic_expansion=approximation,error=error,scaled_t5_error=error*t**5))

out=dict(python=sys.version,scope='scalar finite checks, not full kernel certification',ray_threshold=(3*math.sqrt(3)/(8*pi))**.5,script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),integration=dict(method='Simpson after s=t/(1-u), u in [0,1]',subintervals=8192),formula_rows=rows,derivative_samples=derivative_rows,interior_counterexample=dict(sigma='(delta_(1/2,0)+delta_(-1/2,0))/2',t=counter_t,point=[.5,0],value=interior,exact_interior_upper='-8',exact_boundary_lower='-1'),large_t_samples=large)
print(json.dumps(out,indent=2))
