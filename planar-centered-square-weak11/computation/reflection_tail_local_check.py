#!/usr/bin/env python3
"""Exact parameters for N-independent reflected-periodic tail and local types."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys
from math import isqrt
assert os.path.ismount('/Volumes/data')
def ceilq(v):return -(-v.numerator//v.denominator)
a=argparse.ArgumentParser();a.add_argument('inputs',nargs='+',type=Path);args=a.parse_args()
rows=[];base=None;common_catalog=None;common_labels=None
for path in args.inputs:
    blob=path.read_bytes();inp=json.loads(blob);px=F(inp['Px']);py=F(inp['Py']);lam=F(inp['lam'])
    pts=[tuple(F(x) for x in p) for p in inp['pts']]
    assert len(pts)==204 and len({(x,y) for x,y,m in pts})==204
    assert lam==1 and all(0<=x<px and 0<=y<py and m>0 for x,y,m in pts)
    M=sum((m for x,y,m in pts),F(0));rho=M/(px*py)
    if base is None:base=(px,py,lam,M)
    else:assert base==(px,py,lam,M)
    R0=F(10);Rstar=F(15);eta=F(1,100)
    gold=rho*(1+px/(2*R0))*(1+py/(2*R0))
    gref=rho*(1+px/Rstar)*(1+py/(2*Rstar))
    gsm=rho*(1+(px+eta)/Rstar)*(1+(py/2+eta)/Rstar)
    assert gold<1 and gref<1 and gsm<1
    assert gold==F(23820787253964867284592995960140706111190198205884874124835164741,24052915052186143128797603913750000000000000000000000000000000000)
    cols=ceilq(R0/px);rowspan=ceilq(2*R0/py);rows_count=rowspan+1
    qx=ceilq(Rstar/px)+1;qy=ceilq(Rstar/py)+1
    stable=max(cols,rows_count);poly=max(stable,2*qx+1,2*qy+1)
    assert (cols,rows_count,qx,qy,stable,poly)==(2,5,4,4,5,9)
    xmax=max(x for x,y,m in pts);ymax=max(y for x,y,m in pts)
    xmincut=xmax-R0
    eligible=sum(1 for i in range(-cols+1,1) for x,y,m in pts if x+i*px>xmincut)*rows_count
    assert (qx-1)*px>=Rstar and qx*px>Rstar and qy*py>Rstar
    # Radius/mass bootstrap for every competing right endpoint u>xmax.
    boot=[];owner_R=R0
    for step in range(20):
        nc=ceilq(owner_R/px);nr=ceilq(2*owner_R/py)+1
        row_mass=sum((m for i in range(1-nc,1) for x,y,m in pts if x+i*px>xmax-owner_R),F(0))
        bound=nr*row_mass;v=bound*1000**2/4
        upper=F(isqrt(v.numerator//v.denominator)+1,1000)
        assert 4*upper*upper>bound
        boot.append(dict(R=str(owner_R),columns=nc,rows=nr,mass_bound=str(bound),strict_radius_upper=str(upper),positive_square_margin=str(4*upper*upper-bound)))
        if upper>=owner_R:break
        owner_R=upper
    else:raise AssertionError('bootstrap did not stabilize')
    assert owner_R==F(2051,1000) and nc==1 and nr==2
    labels=[i for i,(x,y,m) in enumerate(pts) if x>xmax-owner_R]
    catalog=[[str(x),str(y+j*py),str(m)] for j in range(nr) for i,(x,y,m) in enumerate(pts) if i in labels]
    catalog_labels=[[i,j] for j in range(nr) for i in labels]
    assert len(labels)==72 and len(catalog)==144
    if common_catalog is None:common_catalog=catalog;common_labels=catalog_labels
    else:assert common_catalog==catalog and common_labels==catalog_labels
    rows.append(dict(input=str(path.resolve()),sha256=hashlib.sha256(blob).hexdigest(),bytes=len(blob),atoms_per_cell=len(pts),Px=str(px),Py=str(py),mass=str(M),rho=str(rho),rho_decimal=float(rho),xmax=str(xmax),ymax=str(ymax),R0=str(R0),old_tail=str(gold),old_tail_decimal=float(gold),Rstar=str(Rstar),reflected_tail=str(gref),reflected_tail_decimal=float(gref),reflected_tail_margin=str(1-gref),cloud_eta_upper=str(eta),cloud_reflected_tail=str(gsm),cloud_tail_decimal=float(gsm),axis_catalog_coarse_columns=cols,axis_catalog_coarse_rows=rows_count,axis_catalog_coarse_filtered_atoms=eligible,owner_radius_bound=str(owner_R),owner_radius_bootstrap=boot,axis_catalog_columns=nc,axis_catalog_rows=nr,axis_catalog_full_atoms=nc*nr*len(pts),axis_catalog_filtered_atoms=len(catalog),axis_catalog_retained_labels=labels,axis_stable_from_N=max(nc,nr),receiver_qx=qx,receiver_qy=qy,polynomial_stable_from_N=poly))
s=Path(__file__).resolve();out=s.parent/'results/reflection-tail-local/exact-parameters.json';out.parent.mkdir(parents=True,exist_ok=True)
catalog_path=s.parent/'inputs/reflection_axis_common_144.json';catalog_path.parent.mkdir(exist_ok=True)
catalog_path.write_text(json.dumps(dict(kind='finite_axis_catalog',lam='1',source_inputs=[r['sha256'] for r in rows],owner_radius_bound='2051/1000',canonical_xmax=rows[0]['xmax'],pts=common_catalog,source_label_and_row=common_labels),ensure_ascii=False,indent=2)+'\n')
catalog_blob=catalog_path.read_bytes()
result=dict(status='PASS',python=sys.version,arithmetic='Fraction exact input and inequalities',random_seed=None,floating_tolerance=None,scope='Tail/local parameter checks only; no N-by-N source construction or halo enumeration',script_sha256=hashlib.sha256(s.read_bytes()).hexdigest(),inputs=rows,common_axis_catalog=dict(path=str(catalog_path),bytes=len(catalog_blob),sha256=hashlib.sha256(catalog_blob).hexdigest(),atoms=len(common_catalog),both_inputs_identical=True))
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='PASS',output=str(out),cases=[dict(input=r['input'],rho=r['rho_decimal'],old_tail=r['old_tail_decimal'],reflected_tail=r['reflected_tail_decimal'],cloud_tail=r['cloud_tail_decimal'],catalog=[r['axis_catalog_columns'],r['axis_catalog_rows']],filtered_atoms=r['axis_catalog_filtered_atoms'],axis_N=r['axis_stable_from_N'],polynomial_N=r['polynomial_stable_from_N']) for r in rows]),ensure_ascii=False,indent=2))
