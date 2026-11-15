#!/usr/bin/env python3
"""Exact R59 direct collar coefficient with fixed total strip mass."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def iv(z):return {'lower':str(z[0]),'upper':str(z[1]),'lower_float_display':float(z[0]),'upper_float_display':float(z[1])}

def main():
    ap=argparse.ArgumentParser()
    for k in ['input','left','right','prior','output']:ap.add_argument('--'+k,required=True)
    a=ap.parse_args();assert os.path.ismount('/Volumes/data')
    paths={k:Path(getattr(a,k)) for k in ['input','left','right','prior']};data={k:json.loads(p.read_text()) for k,p in paths.items()}
    ih=sha(paths['input']);left,right=data['left'],data['right'];prior=data['prior']
    assert prior['input_sha256']==ih
    assert left['status']==right['status']=='PASS' and left['mode']==right['mode']=='area'
    assert left['side']=='left' and right['side']=='right' and F(left['t'])==F(right['t'])
    assert left['input_sha256']==right['input_sha256']==ih
    assert left['script_sha256']==right['script_sha256']
    t=F(left['t']);assert 0<=t<=1
    source=data['input'];m=sum(F(p[2]) for p in source['pts']);assert F(prior['mass_per_cell'])==m
    px=F(source['Px']);X=max(F(p[0]) for p in source['pts'])
    assert F(left['anchor'])==X-3*px and F(right['anchor'])==-3*px
    assert F(left['tail_bound'])<1 and F(right['tail_bound'])<1
    Arec=prior['coefficient_intervals']['A'];A=F(Arec['lower']),F(Arec['upper'])
    L=F(left['true_axis_area_lower']),F(left['true_axis_area_upper']);R=F(right['true_axis_area_lower']),F(right['true_axis_area_upper'])
    BL=L[0]-3*A[1],L[1]-3*A[0];BR=R[0]-3*A[1],R[1]-3*A[0]
    K=BL[0]+BR[0],BL[1]+BR[1]
    J7=(A[0]+L[0]+R[0])/(7*m),(A[1]+L[1]+R[1])/(7*m)
    sign='strict_positive' if K[0]>0 else 'strict_negative' if K[1]<0 else 'undetermined'
    states={row['state']:row['area'] for row in prior['sources']}
    gr=prior['coefficient_intervals']['gamma_rational_axis'];er=F(prior['axis_error_source']['gamma_error'])
    BR0=(F(states['right']['lower'])-2*A[1]+F(gr['lower'])-er,F(states['right']['upper'])-2*A[0]+F(gr['upper'])+er)
    selected_mass=sum(F(p[2]) for p in source['pts'] if F(p[0])==0);added=(1-t)*selected_mass
    costgain=(BR[0]-BR0[1]-A[1]*added/m,BR[1]-BR0[0]-A[0]*added/m)
    if t==F(1,2):assert costgain[0]>F(1,10) and K[1]<-F(14,5)
    result={'status':'PASS','python':sys.version,'script_sha256':sha(__file__),'input_sha256':ih,'t':str(t),'u':str(1-t),'mass_change':'0','cell_mass':str(m),'boundary_program_sha256':left['script_sha256'],'intervals':{k:iv(z) for k,z in [('A',A),('Lambda_L',L),('Lambda_R',R),('B_L',BL),('B_R',BR),('K',K),('J_N7',J7)]},'coefficient_sign':sign,'right_local_mass_costed_increment':iv(costgain),'right_added_mass':str(added),'right_increment_formula':'BR(u)-BR(0)-(A/m)*u*wI; local objective gain only, not a new global positiveK','formula':'J_N=Jper+K/(mN), all integerN>=7, this fixed t and u=1-t','scope':'One fully certified endpoint allocation, not optimization or exclusion of whole continuous t-family; threshold1 only.','random_seed':None,'floating_tolerance':None,'sources':[{'role':k,'path':str(p),'sha256':sha(p),'size':p.stat().st_size} for k,p in paths.items()]}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'PASS','t':str(t),'sign':sign,'K':iv(K),'J_N7':iv(J7)},indent=2))
if __name__=='__main__':main()
