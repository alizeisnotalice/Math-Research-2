"""Gaussian transfer near the critical square-root oscillation scale."""
import argparse,json,os,platform,time
from pathlib import Path
import mpmath as mp


def solve(k,sign):
    with mp.workdps(90):
        d=mp.mpf(2)
        radius=mp.pi*(k+(0 if sign==1 else mp.mpf('.5')))
        t=radius**2
        assert 1-t<-20
        peak=mp.sqrt(d/(2*mp.pi))
        def y_transform(y):
            h=y/(2*radius)
            v=mp.sqrt(1+h*h)
            delta=y+y*y/(2*radius*(v+1))
            u=v+h
            return peak*mp.exp(-d*y*y/2)*(1+sign*mp.cos(delta)/2)/(u**3*v)
        def shift(y):
            h=y/(2*radius)
            v=mp.sqrt(1+h*h)
            return y/2+y*y/(4*radius*(v+1))
        def root_shift(s):
            r=radius+s
            y=s*(2*radius+s)/r
            return 2*peak*mp.exp(-d*y*y/2)*(1+sign*mp.cos(2*s)/2)*(radius/r)**4
        a=[-20,-14,-8,-5,-3,-1,0,1,3,5,8,14,20]
        b=[-18,-12,-7,-4,-2,-1,0,1,2,4,7,12,18]
        result=mp.quad(y_transform,a,method='tanh-sinh',maxdegree=9)
        check=mp.quad(root_shift,[shift(y) for y in b],method='gauss-legendre',maxdegree=8)
        pred=1+sign*mp.exp(-1/(2*d))/2
        tail=mp.mpf(27)/4*401*mp.exp(-400)
        check_tail=mp.mpf(27)/4*325*mp.exp(-324)
        difference=abs(result-check)
        assert difference<mp.mpf('1e-80')
        var=1/d
        co=var+sign*mp.exp(-var/2)*(var-var**2+mp.mpf(3)/8*(3*var**2-var**3)
                    -(3*var**2-6*var**3+var**4)/32)/2
        row=dict(k=k,phase_sign=sign,dps=mp.mp.dps,d=d,t=t,radius=radius,result=result,check=check,
                 prediction=pred,finite_minus_prediction=result-pred,difference_between_methods=difference,
                 direct_tail_upper=tail,alternate_tail_upper=check_tail,
                 t_scaled_error=t*(result-pred),predicted_first_finite_correction=co,
                 error_after_finite_correction=result-pred-co/t)
        return {key:mp.nstr(value,62) if isinstance(value,mp.mpf) else value for key,value in row.items()}


def main():
    p=argparse.ArgumentParser();p.add_argument('k',type=int);p.add_argument('sign',type=int,choices=(-1,1));p.add_argument('out',type=Path);a=p.parse_args()
    begun=time.time();cpu=time.process_time()
    rec=solve(a.k,a.sign)
    rec.update(start_epoch=begun,end_epoch=time.time(),cpu_seconds=time.process_time()-cpu,pid=os.getpid(),python=platform.python_version(),mpmath=mp.__version__)
    a.out.write_text(json.dumps(rec,indent=2)+'\n')

if __name__=='__main__':main()
