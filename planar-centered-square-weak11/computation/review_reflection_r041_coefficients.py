#!/usr/bin/env python3
"""Independent Fraction aggregation and independent cheap axis-catalog count.

Reads actual ops state JSON; never computes a real-template halo area.
"""
from fractions import Fraction as F
from math import ceil
from pathlib import Path
import hashlib,json,os,platform,sys


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilq(q):return -(-q.numerator//q.denominator)


def independent_catalog(pts,px,py,axis,state):
    new=state.endswith('-new');R=F(15 if new else 10)
    end=state not in ['seam-old','seam-new']
    y0,y1=(-3*py,3*py) if end else (F(0),py)
    points=set();tags=set()
    for mirrored in [False,True] if new else [False]:
        # All original source x<=xmax<axis: only these lower preimage
        # bounds matter. Unlike the author's generator we do not derive
        # both physical-window bounds or clip an initially finite i range.
        xlower=axis-R if mirrored else axis-15-R
        for t,(alpha,beta,w) in enumerate(pts):
            imin=ceilq((xlower-alpha)/px)
            jmin=ceilq((y0-R-beta)/py);jmax=(y1+R-beta)//py
            if state.startswith('seam-bottom'):jmin=max(jmin,0)
            if state.startswith('seam-top'):jmax=min(jmax,-1)
            for i in range(imin,1):
                original_x=alpha+i*px
                x=2*axis-original_x if mirrored else original_x
                for j in range(jmin,jmax+1):
                    points.add((x,beta+j*py));tags.add((mirrored,t,i,j))
    nx=len({p[0] for p in points});ny=len({p[1] for p in points})
    return len(points),nx,ny,tags


def main():
    assert os.path.ismount('/Volumes/data')
    root=Path('/Volumes/data/CodexProjects/T2/.local/workspaces')
    ops=root/'pcw-operations-terra/reflection'
    author=root/'pcw-astra-reflection/planar-centered-square-weak11'
    ip=root/'pcw-vertical-compute/planar-centered-square-weak11/computation/inputs/r041_final_0036.json'
    sp=ops/'r041-seam.json';ap=author/'computation/results/reflection-edge/r041-axis-error.json'
    data=json.loads(ip.read_text());seam=json.loads(sp.read_text());axis_cert=json.loads(ap.read_text())
    input_hash=sha(ip);seam_hash=sha(sp)
    assert seam['input_sha256']==axis_cert['input_sha256']==input_hash
    assert axis_cert['seam_sha256']==seam_hash
    lo,hi=F(seam['c_edge_lower']),F(seam['c_edge_upper'])
    assert F(seam['rational_reflection_axis'])==hi
    delta=hi-lo;assert 0<delta<=F(1,1<<64)
    px,py=F(data['Px']),F(data['Py']);pts=[tuple(map(F,p)) for p in data['pts'] if F(p[2])>0]
    assert lo>max(p[0] for p in pts) and len(pts)==204
    mass=sum(p[2] for p in pts)
    edge_hash=sha(author/'computation/reflection_edge_states.py')
    corner_hash=sha(author/'computation/reflection_edge_corners.py')
    union_hash=sha(author/'computation/periodic_halo.py')
    state_names=['bulk','left','right','bottom','top','seam-old','seam-new',
       'corner-lb','corner-lt','corner-rb','corner-rt',
       'seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new']
    states={};bindings=[]
    for s in state_names:
        p=ops/f'r041-{s}.json';r=json.loads(p.read_text())
        assert r['status']=='PASS' and r['state']==s and r['input_sha256']==input_hash
        iscorner=s.startswith('corner-') or s.startswith('seam-bottom-') or s.startswith('seam-top-')
        assert r['script_sha256']==(corner_hash if iscorner else edge_hash)
        assert r['union_module_sha256']==union_hash
        assert F(r['cell_mass'])==mass and F(r['lam'])==1
        assert F(r['inner_area'])<=F(r['outer_area']) and F(r['tail_bound'])<1
        if s.startswith('seam'):
            assert F(r['rational_axis'])==hi and r['seam_sha256']==seam_hash
        states[s]=(F(r['inner_area']),F(r['outer_area']))
        bindings.append({'state':s,'path':str(p),'sha256':sha(p)})
    A=states['bulk']
    E=tuple(sum(states[s][i] for s in ['left','right','bottom','top']) for i in [0,1])
    C=tuple(sum(states[s][i] for s in ['corner-lb','corner-lt','corner-rb','corner-rt']) for i in [0,1])
    gamma=(states['seam-new'][0]-states['seam-old'][1],states['seam-new'][1]-states['seam-old'][0])
    Gamma=tuple(sum(states['seam-'+s+'-new'][i]-states['seam-'+s+'-old'][1-i] for s in ['bottom','top']) for i in [0,1])
    L=(E[0]-8*A[1]+gamma[0],E[1]-8*A[0]+gamma[1])
    F6=(2*E[0]-32*A[1]+C[0]+Gamma[0],2*E[1]-32*A[0]+C[1]+Gamma[1])
    # Independently regenerate lower-axis common catalog counts.
    own_rows=[];stored={r['state']:r for r in axis_cert['states']};err={}
    for s in ['seam-old','seam-new','seam-bottom-old','seam-bottom-new','seam-top-old','seam-top-new']:
        atoms,nx,ny,tags_lo=independent_catalog(pts,px,py,lo,s)
        _,_,_,tags_hi=independent_catalog(pts,px,py,hi,s)
        assert tags_hi<=tags_lo
        K=nx*(nx+1)//2*ny*(ny+1)//2
        end=s not in ['seam-old','seam-new'];R=F(15 if s.endswith('-new') else 10)
        height=min(6*py if end else py,2*R)
        e=2*delta*height*K;rr=stored[s]
        assert (atoms,nx,ny,K)==(rr['atoms'],rr['nx'],rr['ny'],rr['all_bbox_upper_count'])
        assert height==F(rr['height_bound']) and e==F(rr['area_change_upper'])
        err[s]=e;own_rows.append({'state':s,'atoms':atoms,'nx':nx,'ny':ny,'K':K,'height':str(height),'error':str(e)})
    eg=err['seam-new']+err['seam-old']
    eG=sum(e for s,e in err.items() if s not in ['seam-new','seam-old'])
    assert eg==F(axis_cert['gamma_change_upper']) and eG==F(axis_cert['Gamma_change_upper'])
    Ltrue=(L[0]-eg,L[1]+eg);F6true=(F6[0]-eG,F6[1]+eG)
    assert Ltrue[1]<-9 and F6true[1]<-53
    # Independent label-only bound, with no data-coordinate catalog count.
    cols=ceil(30/px)+1;rows=ceil(30/py)+2;endrows=ceil(30/py)+7;n=len(pts)
    def boxes(x,y):return x*(x+1)//2*y*(y+1)//2
    coarse_kg=boxes(2*n*cols,n*rows)+boxes(n*cols,n*rows)
    coarse_ke=2*(boxes(2*n*cols,n*endrows)+boxes(n*cols,n*endrows))
    coarse_g=2*delta*py*coarse_kg;coarse_G=12*delta*py*coarse_ke
    assert coarse_g<F(1,100000) and coarse_G<F(1,1000)
    assert L[1]+coarse_g<-9 and F6[1]+coarse_G<-53
    def encode(v):return list(map(str,v))
    output={'status':'PASS','scope':'Independent exact aggregation and cheap catalog counts; source area JSON reused, not recomputed.',
       'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'single_process':True,
       'arithmetic':'Fraction','random_seed':None,'floating_tolerance':None,
       'script_sha256':sha(Path(__file__)),'input_sha256':input_hash,'seam_sha256':seam_hash,
       'axis_error_json_sha256':sha(ap),'edge_kernel_sha256':edge_hash,'corner_kernel_sha256':corner_hash,
       'union_module_sha256':union_hash,'bindings':bindings,'delta_upper':str(delta),
       'A_interval':encode(A),'E_interval':encode(E),'Ccorner_interval':encode(C),
       'gamma_rational_axis_interval':encode(gamma),'Gamma_rational_axis_interval':encode(Gamma),
       'L_rational_axis_interval':encode(L),'F6_rational_axis_interval':encode(F6),
       'independent_catalog_rows':own_rows,'gamma_axis_error':str(eg),'Gamma_axis_error':str(eG),
       'L_true_axis_interval':encode(Ltrue),'F6_true_axis_interval':encode(F6true),
       'label_only_bound':{'column_count':cols,'bulk_rows':rows,'end_rows':endrows,
          'gamma_K_upper':coarse_kg,'Gamma_K_upper':coarse_ke,'gamma_error':str(coarse_g),'Gamma_error':str(coarse_G)},
       'clean_rational_sign_certificates':{'L_true_axis_upper':'<-9','F6_true_axis_upper':'<-53'},
       'conclusion':'For the specified fixed-phase R41 template and every N>=6, both rational outward-axis reflection and true-axis reflection have ratio strictly below the exact periodic bulk A/m. N=1..5 and other templates/phases remain open.'}
    out=Path(__file__).parent/'results/reflection-edge-state-review/r041-independent-coefficients.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','state_json_count':len(bindings),'L_true_upper_float':float(Ltrue[1]),
       'F6_true_upper_float':float(F6true[1]),'axis_gamma_error':float(eg),'axis_Gamma_error':float(eG),
       'coarse_gamma_error':float(coarse_g),'coarse_Gamma_error':float(coarse_G),'output':str(out)},ensure_ascii=False))

if __name__=='__main__':main()
