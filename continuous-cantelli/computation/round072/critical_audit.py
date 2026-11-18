"""Independent closed-Gaussian interval audit of the Cantor renewal limit."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=63
omega=2*mp.pi

def kernel(a,T):
    z=mp.j*omega/(2*mp.sqrt(a))
    return mp.sqrt(mp.pi)/(2*mp.sqrt(a))*mp.exp(-omega**2/(4*a))*(mp.erf(mp.sqrt(a)*T-z)-mp.erf(-z))

def renewal(b,a,k=None,terms=145):
    H=mp.mpf(b-2)/2
    phi=-mp.mpf(1)
    scale=mp.sqrt(mp.pi)/(2*mp.sqrt(a))
    weighted=[]
    for l in range(1,terms+1):
        t=H/b**l
        factor=mp.mpf(2)**l
        weighted.append((factor*scale*mp.erf(mp.sqrt(a)*t),
                         factor*phi*mp.re(mp.exp(-mp.j*omega*t)*kernel(a,t))))
        phi*=mp.cos(mp.pi*(b-1)/mp.mpf(b)**l)
    for s in range(terms if k is None else k):
        t=H*b**s
        q=mp.mpf(2)**(-s)
        weighted.append((q*scale*mp.erf(mp.sqrt(a)*t),q*mp.re(kernel(a,t))))
    den=mp.fsum(x[0] for x in weighted)
    num=mp.fsum(x[1] for x in weighted)
    return den,num,num/den

def main():
    out=[]
    for b in (3,4,6,7):
        for v in ('.003','.03','.3','1','8','64'):
            a=mp.mpf(v)
            de,nu,value=renewal(b,a)
            surrogate=mp.re(mp.exp(-mp.j*mp.pi)*mp.fprod(mp.cos(mp.pi*(b-1)/mp.mpf(b)**j)
                       for j in range(1,145)))*mp.exp(-mp.pi**2/a)
            out.append(dict(base=b,a=v,critical_ratio=str(value),unnormalized_scaled=str(nu),
                            Gaussian_convolution_surrogate=str(surrogate)))
    lo=mp.mpf('.003');hi=mp.mpf(64)
    assert renewal(4,lo)[1]<0<renewal(4,hi)[1]
    for _ in range(125):
        mid=(lo+hi)/2
        if renewal(4,mid)[1]>0: hi=mid
        else: lo=mid
    a=(lo+hi)/2
    finite=[]
    tail_coefficient=mp.sqrt(mp.pi)/mp.sqrt(a)*mp.exp(-mp.pi**2/a)
    for k in (4,8,12,16,22):
        den,g,val=renewal(4,a,k=k)
        finite.append(dict(k=k,normalized_gibbs=str(val),missing_tail_scaled=str(-2**k*g),
                           predicted=str(tail_coefficient)))
    result=dict(digits=mp.mp.dps,series_terms=145,cases=out,illustrative_base4_root=str(a),root_bracket=[str(lo),str(hi)],
                finite_root=finite)
    Path('critical_audit.json').write_text(json.dumps(result,indent=2)+"\n")
    print("b4 root",a)
    print("finite",finite[-1])

if __name__=='__main__': main()
