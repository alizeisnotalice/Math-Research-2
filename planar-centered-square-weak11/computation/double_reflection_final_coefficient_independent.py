#!/usr/bin/env python3
"""Independent R54 final Fraction aggregation; uses completed area certificates."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json,os,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    assert os.path.ismount('/Volumes/data')
    base=Path(__file__).resolve().parents[1]
    work=Path('/Volumes/data/CodexProjects/T2/.local/workspaces')
    ops=work/'pcw-operations-terra'
    core=work/'pcw-astra-reflection/planar-centered-square-weak11'
    paths={name:ops/'reflection'/('r041-'+name+'.json') for name in ('bulk','left','right')}
    paths.update({name:ops/'r054'/('r041-'+name+'.json') for name in ('left-axis','leftseam-old','leftseam-new')})
    paths['right-certificate']=base/'computation/results/reflection-edge-state-review/r041-independent-coefficients.json'
    paths['left-axis-error']=core/'computation/results/double-reflection/r041-left-axis-error.json'
    paths['transform']=core/'computation/results/double-reflection/r041-transform.json'
    paths['transformed-input']=core/'computation/inputs/r054_r041_left_reflected.json'
    paths['original-input']=work/'pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json'
    d={k:read(v) for k,v in paths.items()}
    original=d['original-input'];flipped=d['transformed-input'];tr=d['transform']
    X=max(F(p[0]) for p in original['pts']);px=F(original['Px'])
    assert X==F(tr['X']) and 0<=X<px
    assert sorted((X-F(x),F(y),F(w)) for x,y,w in original['pts'])==sorted(tuple(map(F,p)) for p in flipped['pts'])
    assert sha(paths['transformed-input'])==d['left-axis']['input_sha256']==d['leftseam-old']['input_sha256']==d['leftseam-new']['input_sha256']
    assert sha(paths['left-axis'])==d['leftseam-old']['seam_sha256']==d['leftseam-new']['seam_sha256']==d['left-axis-error']['seam_sha256']
    delta=F(d['left-axis']['c_edge_interval_width']);assert delta==F(d['left-axis-error']['delta_upper'])
    states={s['state']:s for s in d['left-axis-error']['states']}
    errors=[]
    for name in ('seam-old','seam-new'):
        s=states[name];boxes=comb(s['nx']+1,2)*comb(s['ny']+1,2)
        assert boxes==s['all_bbox_upper_count']
        error=2*delta*F(s['height_bound'])*boxes
        assert error==F(s['area_change_upper']);errors.append(error)
    eg=sum(errors,F(0));assert eg==F(d['left-axis-error']['gamma_change_upper'])
    old,new=d['leftseam-old'],d['leftseam-new']
    gl=(F(new['inner_area'])-F(old['outer_area'])-eg,F(new['outer_area'])-F(old['inner_area'])+eg)
    left,right,bulk=d['left'],d['right'],d['bulk']
    bx=(F(left['inner_area'])+F(right['inner_area'])-4*F(bulk['outer_area']),F(left['outer_area'])+F(right['outer_area'])-4*F(bulk['inner_area']))
    rc=d['right-certificate'];gr0=tuple(map(F,rc['gamma_rational_axis_interval']));er=F(rc['gamma_axis_error'])
    gr=(gr0[0]-er,gr0[1]+er)
    K=(bx[0]+gl[0]+gr[0],bx[1]+gl[1]+gr[1])
    assert K[1]<-F(12,5)
    out={'status':'PASS','scope':'Exact aggregation of completed certificates, transform and axis-error arithmetic; no rerun of full bbox area kernels.',
         'python':sys.version,'arithmetic':'Fraction','random_seed':None,'floating_tolerance':None,
         'gammaL_true_interval':list(map(str,gl)),'gammaR_true_interval':list(map(str,gr)),
         'bx_interval':list(map(str,bx)),'K_interval':list(map(str,K)),
         'display_gammaL':list(map(float,gl)),'display_K':list(map(float,K)),
         'exact_upper_bound':'-12/5','exact_negative_margin':str(-F(12,5)-K[1]),
         'consequence':'For every N>=7, J_N(lambda=1)<A/m-12/(5mN); for fixed N the specified finite-L1 sequence eventually has J(f_j;1)<A/m-6/(5mN). No sup-lambda bound.',
         'source_hashes':{str(p):sha(p) for p in paths.values()},'script_sha256':sha(Path(__file__))}
    target=base/'computation/results/double-reflection-closure/final-coefficient-independent.json'
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','gammaL':list(map(float,gl)),'K':list(map(float,K)),
                     'K_upper_lt_minus_12_over_5':True,'output':str(target),'sha256':sha(target)}))
if __name__=='__main__':main()
