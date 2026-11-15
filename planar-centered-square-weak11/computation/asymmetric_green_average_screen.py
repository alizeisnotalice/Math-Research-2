#!/usr/bin/env python3
"""Deterministic scalar screen; samples are necessary constraints, not D certificates."""
import os,sys,math,cmath,json,time,hashlib,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def green(z,s,b):
    u=z/s
    if b==0:w=u
    else:w=2*u/(1+cmath.sqrt(1+4*b*u))
    r=abs(w)
    return max(0.,-math.log(r))/(2*math.pi)

def avg(z,s,b):return sum(green(z*(1j**k),s,b) for k in range(4))/4

def evaluate(s,b):
    vals=[avg(complex(1-t,t),s,b) for t in ANGLES]
    m=min(vals)
    return (math.inf if m<=0 else s*s*(1+b*b)/(8*m)),ANGLES[vals.index(m)],m

ANGLES=[j/64 for j in range(33)]
if __name__=='__main__':
    assert os.path.ismount('/Volumes/data')
    start=time.monotonic();best=[];byb=[];count=0
    for j in range(32):
        b=j/64;rows=[]
        for k in range(48,225):
            s=k/64;cost,t,m=evaluate(s,b);count+=1
            rows.append({'s':s,'b':b,'sample_necessary_cost':cost,'argmin_sample_t':t,'min_G_sample':m})
        byb.append(min(rows,key=lambda r:r['sample_necessary_cost']))
        best+=sorted(rows,key=lambda r:r['sample_necessary_cost'])[:2]
    best=sorted(best,key=lambda r:r['sample_necessary_cost'])[:12]
    # Deterministic local refinement near best coarse s, including every sampled b.
    refined=[]
    for row in byb:
        s0=row['s'];b=row['b'];local=[]
        for j in range(-32,33):
            s=s0+j/4096;v,t,m=evaluate(s,b);count+=1
            local.append({'s':s,'b':b,'sample_necessary_cost':v,'argmin_sample_t':t,'min_G_sample':m})
        refined.append(min(local,key=lambda r:r['sample_necessary_cost']))
    out={'status':'COMPLETED_numerical_screen_only','python':sys.version,'platform':platform.platform(),'workers':1,'elapsed_seconds':time.monotonic()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'family':'D4 average of Green kernels of F(w)=s(w+bw^2), b in [0,1/2), group rotations 1,i,-1,-i','mass_formula':'integral K_unscaled=s^2*(1+b^2)/4','parameters':{'coarse_b':'j/64, j=0..31','coarse_s':'k/64,k=48..224','boundary_t':'j/64,j=0..32; z=1-t+it','refinement':'s0+j/4096,j=-32..32, each coarse b','evaluations':count,'random_seed':None},'top_coarse':best,'refined_best_by_b':refined,'best_refined':min(refined,key=lambda r:r['sample_necessary_cost']),'circle_exact_optimum':math.pi*math.e/2,'interpretation':'Grid sample min G >= true boundary min G, hence printed ratios are sampled necessary costs, not feasible upper bounds; no continuous parameter optimality or boundary positivity certificate.'}
    dest=ROOT/'computation/results/asymmetric-green-average-screen.json';dest.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'evaluations':count,'elapsed_seconds':out['elapsed_seconds'],'best_refined':out['best_refined'],'best_b_positive':min([r for r in refined if r['b']>0],key=lambda r:r['sample_necessary_cost']),'output':str(dest)},indent=2))
