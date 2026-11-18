"""Bivariate real-coefficient Taylor models with a holomorphic polydisc remainder."""
from math import factorial

from flint import arb


def powers(x,d):
    row=[arb(1)]
    for _ in range(d):
        row.append(row[-1]*x)
    return row


class Model:
    degree=8

    def __init__(self,value=0,error=0):
        if isinstance(value,dict):
            self.poly=value
        else:
            v=arb(value)
            self.poly={(0,0):v} if v else {}
        self.error=arb(error).abs_upper()
        self.bound=sum((x.abs_upper() for x in self.poly.values()),arb(0)).abs_upper()

    @staticmethod
    def promote(a):
        return a if isinstance(a,Model) else Model(a)

    def __add__(self,a):
        a=self.promote(a)
        out=self.poly.copy()
        for k,v in a.poly.items():
            out[k]=out.get(k,0)+v
        return Model(out,self.error+a.error)

    __radd__=__add__

    def __neg__(self):
        return Model({k:-v for k,v in self.poly.items()},self.error)

    def __sub__(self,a):
        return self+-self.promote(a)

    def __rsub__(self,a):
        return self.promote(a)+-self

    def __mul__(self,a):
        if not isinstance(a,Model):
            a=arb(a)
            return Model({k:v*a for k,v in self.poly.items()},self.error*a.abs_upper())
        polynomial={}
        omitted=arb(0)
        for (i,j),v in self.poly.items():
            for (k,l),w in a.poly.items():
                if i+j+k+l<=self.degree:
                    loc=(i+k,j+l)
                    polynomial[loc]=polynomial.get(loc,0)+v*w
                else:
                    omitted+=v.abs_upper()*w.abs_upper()
        return Model(polynomial,omitted+self.error*a.bound+a.error*self.bound+self.error*a.error)

    __rmul__=__mul__

    def __truediv__(self,a):
        return self*a.inv() if isinstance(a,Model) else self*(1/arb(a))

    def __rtruediv__(self,a):
        return self.inv()*a

    def __pow__(self,k):
        assert isinstance(k,int) and k>=0
        p=Model(1)
        for _ in range(k):
            p=p*self
        return p

    def separate(self):
        c=self.poly.get((0,0),arb(0))
        rest=Model({key:v for key,v in self.poly.items() if key!=(0,0)},self.error)
        return c,rest,(rest.bound+rest.error).abs_upper()

    def compose(self,rest,coeff,tail):
        result=Model(coeff[-1])
        for c in reversed(coeff[:-1]):
            result=result*rest+c
        return Model(result.poly,result.error+tail)

    def inv(self):
        c,w,r=self.separate()
        lo=c.abs_lower()
        assert r<lo,(r,lo)
        d=self.degree
        rho=r/lo
        w=w*(1/c)
        tail=rho**(d+1)/(lo*(1-rho))
        return self.compose(w,[((-1)**j)/c for j in range(d+1)],tail)

    def sqrt(self):
        c,w,r=self.separate()
        assert c>r,(c,r)
        root=c.sqrt()
        q=r/c.lower()
        w=w*(1/c)
        coeff=[root]
        for j in range(self.degree):
            coeff.append(coeff[-1]*(arb(1)/2-j)/(j+1))
        tail=root.abs_upper()*q**(self.degree+1)/(1-q)
        return self.compose(w,coeff,tail)

    def atan(self):
        c,w,r=self.separate()
        p=1+c*c
        q=r/p.sqrt().lower()
        assert q<1
        real,imag=c/p,1/p
        u,v=arb(1),arb(0)
        coeff=[c.atan()]
        for k in range(1,self.degree+1):
            u,v=u*real-v*imag,u*imag+v*real
            coeff.append((-1)**(k+1)*v/k)
        tail=q**(self.degree+1)/((self.degree+1)*(1-q))
        return self.compose(w,coeff,tail)

    def _trig(self,want_sine):
        c,w,r=self.separate()
        ss,cc=c.sin(),c.cos()
        derivatives=[ss,cc,-ss,-cc] if want_sine else [cc,-ss,-cc,ss]
        coeff=[derivatives[j%4]/factorial(j) for j in range(self.degree+1)]
        tail=r.exp()*r**(self.degree+1)/factorial(self.degree+1)
        return self.compose(w,coeff,tail)

    def sin(self):
        return self._trig(True)

    def cos(self):
        return self._trig(False)

    def evaluate(self,x,y):
        xp,yp=powers(x,self.degree),powers(y,self.degree)
        p=sum((v*xp[i]*yp[j] for (i,j),v in self.poly.items()),arb(0))
        return arb(p,self.error)

    def second(self,i,j,x,y,scale):
        m,n=(2,0) if i==j==0 else (0,2) if i==j==1 else (1,1)
        value=arb(0)
        xp,yp=powers(x,self.degree),powers(y,self.degree)
        for (a,b),v in self.poly.items():
            if a>=m and b>=n:
                f=(a*(a-1) if m==2 else a if m else 1)*(b*(b-1) if n==2 else b if n else 1)
                value+=v*f*xp[a-m]*yp[b-n]
        margin=1-max(x.abs_upper(),y.abs_upper())
        assert margin>0
        # Independent Cauchy circles in the two normalized complex variables.
        extra=self.error*(2 if i==j else 1)/margin**2
        return arb(value,extra.abs_upper())/scale**2
