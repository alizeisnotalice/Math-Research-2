#!/usr/bin/env python3
"""Directed complete-catalog area bound from rational outer axis to true axis."""
import argparse,hashlib,json,os,sys
from fractions import Fraction as F
from pathlib import Path
import reflection_edge_states as edges
import reflection_edge_corners as corners

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--seam',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data');src=Path(a.input);sp=Path(a.seam);data=json.loads(src.read_text());seam=json.loads(sp.read_text())
    assert hashlib.sha256(src.read_bytes()).hexdigest()==seam['input_sha256']
    px,py,lam=(F(data[k]) for k in ('Px','Py','lam'));pts=[tuple(map(F,p)) for p in data['pts'] if F(p[2])>0]
    lo,hi=F(seam['c_edge_lower']),F(seam['c_edge_upper']);delta=hi-lo
    assert 0<delta<=F(1,1<<64) and lo>max(x for x,y,m in pts)
    rows=[]
    for state in ('seam-old','seam-new','seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new'):
        R=F(15 if state.endswith('-new') else 10)
        iscorner=state not in ('seam-old','seam-new');win=(lo-15,lo,-3*py,3*py) if iscorner else (lo-15,lo,F(0),py)
        module=corners if iscorner else edges
        atoms=module.source_atoms(pts,px,py,win,R,state,lo)
        nx=len({x for x,y in atoms});ny=len({y for x,y in atoms});K=nx*(nx+1)//2*ny*(ny+1)//2
        clipped_height=min(win[3]-win[2],2*R);err=2*delta*clipped_height*K
        rows.append({'state':state,'common_catalog_axis':'lower interval endpoint','atoms':len(atoms),'nx':nx,'ny':ny,'all_bbox_upper_count':K,'height_bound':str(clipped_height),'area_change_upper':str(err),'area_change_upper_float':float(err)})
    by={r['state']:F(r['area_change_upper']) for r in rows}
    gamma=by['seam-new']+by['seam-old'];Gamma=sum(v for k,v in by.items() if k not in ('seam-new','seam-old'))
    result={'status':'PASS','python':sys.version,'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'seam_sha256':hashlib.sha256(sp.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'edge_state_code_sha256':hashlib.sha256(Path(edges.__file__).read_bytes()).hexdigest(),'corner_state_code_sha256':hashlib.sha256(Path(corners.__file__).read_bytes()).hexdigest(),'delta_upper':str(delta),'states':rows,'gamma_change_upper':str(gamma),'Gamma_change_upper':str(Gamma),'gamma_change_upper_float':float(gamma),'Gamma_change_upper_float':float(Gamma),'proof_scope':'In axis-centered coordinates fixed window [-15,0], original and mirror source x move by at most delta. Their order is fixed because axis lower endpoint exceeds all original x. Each bbox receiver endpoint moves at most delta; every nonempty receiver at any axis has radius<R by the uniform tail. Complete catalog at the lower axis is a superset for the entire axis interval. Area symmetric difference per bbox <=2 delta min(window height,2R), including empty/activation/cropping events. Counts include all boxes, not only surviving rectangles.'}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','gamma_change_upper_float','Gamma_change_upper_float')},ensure_ascii=False))
if __name__=='__main__':main()
