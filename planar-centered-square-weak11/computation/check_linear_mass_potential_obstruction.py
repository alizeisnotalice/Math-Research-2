#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,sys,platform

def singleton_area(h,j,l):
    t=l*h
    width=max(F(0),min(h,j*h+t)-max(F(0),j*h-t))
    return 2*t*width

def main():
    assert os.path.ismount('/Volumes/data')
    h=F(1,64); T=F(4); R=F(1); k=65; L=32
    beta={}; tests=0
    for j in range(-k,k+1):
        ratios=[]
        for l in range(1,L+1):
            m=4*l*l*h*h; a=singleton_area(h,j,l)
            ratios.append(a/m); tests+=1
            assert l*h<R and 2*l*h<T
        beta[j]=max(ratios)
        expected=F(1,2*j) if 1<=j<=L else F(1,2*(1-j)) if 1-L<=j<=0 else F(0)
        assert beta[j]==expected
    z=sum(beta.values(),F(0)); H=sum((F(1,l) for l in range(1,L+1)),F(0))
    assert z==H and z>4
    alpha={-k-1:F(0)}
    for j in range(-k,k+1):
        alpha[j]=alpha[j-1]-beta[j]+(z if j==0 else 0)
    assert alpha[k]==0
    checked=0
    for j in range(-k,k+1):
        for l in range(1,L+1):
            m=4*l*l*h*h; a=singleton_area(h,j,l)
            dp=(alpha[j-1]-alpha[j])*m
            assert a-(z*m if j==0 else 0)<=dp
            checked+=1
    levels=[]
    for n in range(2,11):
        hn=F(1,2**n); maxl=2**(n-2); step=F(1,2**(4*n))
        for l in range(1,maxl+1):
            m=4*l*l*hn*hn
            assert m<=F(1,4) and (m/step).denominator==1
            assert l*hn<F(2**n) and 2*l*hn<F(2**n)
        Hn=sum((F(1,l) for l in range(1,maxl+1)),F(0))
        assert Hn>=F(n,2)
        levels.append({'n':n,'test_distance_count':maxl,'harmonic_bound':str(Hn),'simple_bound':str(F(n,2))})
    source=Path(__file__).resolve()
    result={'status':'PASS','python_version':sys.version,'platform':platform.platform(),
      'script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'arithmetic':'fractions.Fraction','seed':None,'tolerance':None,
      'parameters':{'h':str(h),'T':str(T),'R':str(R),'k':k,'mass_count':L},
      'exact_singleton_area_tests':tests,'potential_inequality_tests':checked,'optimal_singleton_only_linear_cost':str(z),
      'lower_than_general_model_scope':False,'scope':'Exact optimum of the singleton-test linear envelope; no multicolumn upper certificate claimed',
      'beta':{str(j):str(b) for j,b in beta.items()},'alpha':{str(j):str(a) for j,a in alpha.items()},'hierarchy_checks':levels}
    out=source.parent/'results/linear-mass-potential/obstruction.json'
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','single_column_area_tests':tests,'potential_tests':checked,'H32':str(z),'H32_gt_4':z>4,'levels':len(levels)},indent=2))
    print(out)
if __name__=='__main__': main()
