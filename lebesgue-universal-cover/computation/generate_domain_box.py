"""Exact-arithmetic outer-domain certificate for the proposed eight test bodies.

It proves only that any placement with hull area < 419/500 lies in the
reported compact parameter box.  It does not bound the area inside that box.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import hashlib, json, math, platform, sys
from pathlib import Path

OUT = Path(__file__).with_name("outer_domain_certificate.json")
TAU = F(419, 500)
W_REGULAR = F(21, 100)
W_FOURIER = F(6, 25)
S_REGULAR = {3:F(7061,10000), 5:F(7064,10000), 7:F(7079,10000), 9:F(7086,10000), 11:F(709,1000), 19:F(887,1250)}
S_FOURIER = F(427,600)  # 283/600 + 6/25
GRID_DEN = 10**15
GRID_N = 512
T = 10**18

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b):
    z=[x*y for x in a for y in b]
    return min(z),max(z)
def scale(a,q): return mul(a,(q,q))
def inv(a):
    assert a[0]>0
    return 1/a[1],1/a[0]
def arctan_inv(q,N):
    v=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F())
    r=F(1,(2*N+1)*q**(2*N+1))
    return (v,v+r) if N%2==0 else (v-r,v)

# Machin's identity gives this certified rational pi interval.
t2=2*F(1,5)/(1-F(1,5)**2); t4=2*t2/(1-t2*t2)
assert t2==F(5,12) and t4==F(120,119) and (t4-F(1,239))/(1+t4*F(1,239))==1
machin=add(scale(arctan_inv(5,32),F(16)),scale(arctan_inv(239,8),F(-4)))
PI=(F(3141592653589793238,T),F(3141592653589793239,T))
assert PI[0]<machin[0]<=machin[1]<PI[1]

@lru_cache(None)
def cos_interval(a):
    mid=(a[0]+a[1])/2
    x=F((mid*T).__floor__(),T)
    rad=max(abs(a[0]-x),abs(a[1]-x))
    assert abs(x)<12
    term=F(1); value=term
    for j in range(1,41):
        term *= -x*x/((2*j-1)*(2*j)); value += term
    rem=abs(x)**82/F(factorial(82))+rad
    return F((value-rem)*T//1,T),F(-((-value-rem)*T//1),T)

def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def hull(P):
    P=sorted(set(P)); ans=[]
    for seq in (P,list(reversed(P))):
        h=[]
        for p in seq:
            while len(h)>1 and cross(h[-2],h[-1],p)<=0: h.pop()
            h.append(p)
        ans += h[:-1]
    return ans
def area(H):
    return F(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(H,H[1:]+H[:1])),2*GRID_DEN*GRID_DEN)

def main():
    # These are merely a deterministic way to choose rational disk points.
    # Every subsequent membership and area statement is exact integer/Fraction arithmetic.
    shrink=1.0-1e-12
    disk=[]
    for j in range(GRID_N):
        disk.append((int(.5*shrink*math.cos(2*math.pi*j/GRID_N)*GRID_DEN), int(.5*shrink*math.sin(2*math.pi*j/GRID_N)*GRID_DEN)))
    assert all(x*x+y*y <= (GRID_DEN//2)**2 for x,y in disk)
    rows=[]
    for n,S in S_REGULAR.items():
        ch=cos_interval(scale(PI,F(1,2*n)))
        cn=cos_interval(scale(PI,F(1,n)))
        Rlo=inv(scale(ch,F(2)))[0]
        g2lo=Rlo*Rlo+W_REGULAR*W_REGULAR+2*Rlo*W_REGULAR*cn[0]
        assert S*S < g2lo
        H=hull(disk+[(int(S*GRID_DEN),0)])
        A=area(H)
        assert A>TAU
        rows.append({'body':f'regular_{n}','translation_coordinate_halfwidth':str(W_REGULAR),'point_radius_lower':str(S),'radius_squared_lower_strict':str(g2lo),'inscribed_polygon_area':str(A),'hull_vertices':len(H)})
    HF=hull(disk+[(int(S_FOURIER*GRID_DEN),0)])
    AF=area(HF)
    assert AF>TAU
    rows.append({'body':'fourier_h=1/2+cos(3t)/80-19cos(5t)/1200','translation_coordinate_halfwidth':str(W_FOURIER),'support_lower':'283/600','point_radius_lower':str(S_FOURIER),'inscribed_polygon_area':str(AF),'hull_vertices':len(HF)})
    result={'status':'verified','scope':'if the eight-body hull has area <419/500, every normalized placement is in this box; no interior branch-and-bound claim','normalization':'D at 0; B3 angle 0; use whole-placement reflection to put rho_B5 in [0,pi/5]','parameter_order':['t_B3','rho_B5','t_B5','rho_B7','t_B7','rho_B9','t_B9','rho_P11','t_P11','rho_P19','t_P19','rho_F','t_F'],'translation_boxes':{'B3':'[-21/100,21/100]^2','B5':'[-21/100,21/100]^2','B7':'[-21/100,21/100]^2','B9':'[-21/100,21/100]^2','P11':'[-21/100,21/100]^2','P19':'[-21/100,21/100]^2','Fourier':'[-6/25,6/25]^2'},'angle_boxes':{'rho_B5':'[0,pi/5]','rho_B7':'[0,2pi/7]','rho_B9':'[0,2pi/9]','rho_P11':'[0,2pi/11]','rho_P19':'[0,2pi/19]','rho_F':'[0,2pi]'},'tau':str(TAU),'disk_polygon_denominator':GRID_DEN,'disk_polygon_points':disk,'regular_checks':rows[:-1],'fourier_check':rows[-1],'pi_interval':[str(PI[0]),str(PI[1])],'python':sys.version,'platform':platform.platform(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'rows':rows,'output':str(OUT)},indent=2))
if __name__=='__main__': main()
