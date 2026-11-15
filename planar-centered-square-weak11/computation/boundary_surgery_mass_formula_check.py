#!/usr/bin/env python3
"""R59 exact mass/accounting prerequisites, not a halo-area t scan."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert os.path.ismount('/Volumes/data')
    src=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json')
    d=json.loads(src.read_text());px=F(d['Px']);py=F(d['Py']);pts=[tuple(map(F,p)) for p in d['pts']]
    I=(0,51);assert all(pts[i][0]==0 for i in I)
    wI=sum((pts[i][2] for i in I),F(0));m=sum(p[2] for p in pts)
    assert 3*px>15 and 6*px>30
    # Minimal N=7 is sufficient for 3-cell boundary windows plus R15
    # source padding, assuming ell<0<r. More cells only enlarge padding.
    assert 7*px>3*px+15 and -7*px+max(p[0] for p in pts)<-3*px-15
    records=[]
    for t in (F(0),F(1,2),F(1)):
        u=1-t
        each=[(7-1+t+u)*pts[i][2] for i in I]
        assert all(mm==7*pts[i][2] for i,mm in zip(I,each))
        total=7*(m-wI)+(7-1+t+u)*wI
        assert total==7*m
        records.append({'t':str(t),'u':str(u),'strip_mass':str(total),'double_cell_mass':str(2*total)})
    out={'status':'PASS','scope':'Actual R41 exposed-label mass conservation and R15 window-padding prerequisites only; no halo t scan.',
         'python':sys.version,'arithmetic':'Fraction','N_min':7,'exposed_labels_zero_based':list(I),
         'exposed_points':[[str(x) for x in pts[i]] for i in I],
         'w_I':str(wI),'m':str(m),'Px':str(px),'Py':str(py),'three_Px_minus_15':str(3*px-15),
         'balanced_cases':records,'general_mass_delta':'(t+u-1)*w_I',
         'general_gain':'[BL(t)+BR(u)-Jper*delta_m]/[N*m+delta_m]',
         'required_axis_conditions':'ell(t)<0<r(u), all positive atoms strictly between axes; uniform full R15 and window catalogs must be certified separately',
         'source_sha256':sha(src),'script_sha256':sha(Path(__file__)),'random_seed':None,'floating_tolerance':None}
    dest=Path(__file__).resolve().parent/'results/boundary-surgery/mass-formula-check.json'
    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','w_I':str(wI),'3Px_minus15':str(3*px-15),'output':str(dest),'sha256':sha(dest)}))
if __name__=='__main__':main()
