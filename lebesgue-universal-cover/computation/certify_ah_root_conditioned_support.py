"""Use the same centering translation through alternative hexagon support pairs at a compatibility root."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_ah_coupled_brackets import support_bound
from certify_ah_projected_jets import geometry


def root_conditioned(ns,j,left,right,angular,begin=0,end=75):
    A=ns['A'];z=A(0);one=A(1);cross=ns['cross'];rotate=ns['rotate'];trig=ns['trig'];support=ns['supports'];
    normals=list(map(ns['av'],ns['d']['n']));assert len(normals)==6
    turn,sign,orient=(300,1,1) if j<34 else (0,1,-1)
    conditions=[]
    for name,u,c in ns['rows']:
        pairs=[]
        for k in range(6):
            n,m=normals[k],normals[(k+1)%6];det=cross(n,m);assert det>0
            alpha,beta=cross(u,m)/det,cross(n,u)/det
            if k==0 or (not alpha<0 and not beta<0):pairs.append((k,alpha,beta))
        assert pairs
        conditions.append((name,u,c,pairs))
    worstroot=worstpose=one;worstrootid=worstposeid=None;changed=0;max_gain=z;witness=None
    for i in range(begin,end):
        a=F(j,100)+F(i,7500);b=a+F(1,7500);E,v,e,bo,_=support_bound(ns,j,a,b)
        for phi,s in ((left,orient),(right,-orient)):
            er=bo([(A(-s*(-1)**k),rotate(normals[0],trig(A(phi+60*k)))) for k in range(3)]);val=(s*ns['F_at'](v,e,phi)-E*er).lower()
            if val<worstroot:worstroot,worstrootid=val,[i,str(phi)]
        for p in range(angular):
            l,r=left+(right-left)*F(p,angular),left+(right-left)*F(p+1,angular)
            mt=trig(A((l+r)/2-turn));sn=trig(ns['box'](l-turn,r-turn));dw=ns['angle_unit']*ns['box']((l-r)/2,(r-l)/2)
            mid=[support(v,e,mt,n) for n in normals];dv=[support(v,e,sn,n,True) for n in normals];un=[rotate(n,sn) for n in normals]
            for name,u,c,pairs in conditions:
                um=support(v,e,mt,u);du=support(v,e,sn,u,True);ur=rotate(u,sn);vals=[]
                for k,alpha,beta in pairs:
                    k1=(k+1)%6
                    base=c-um+alpha*(mid[k]-A('0.5'))+beta*(mid[k1]-A('0.5'))+dw*(-du+alpha*dv[k]+beta*dv[k1]);rough=(base-A('1.07')*(one+abs(alpha)+abs(beta))*E).lower()
                    val=rough if rough>A('.00028') else (base-E*bo([(one,ur),(-alpha,un[k]),(-beta,un[k1])])).lower();vals.append((val,k))
                selected,kbest=max(vals,key=lambda pair:pair[0]);old=vals[0][0]
                if old<=0 and selected>0:changed+=1
                # Compare a negative legacy bound with a certified positive alternative.
                if old<=0 and selected>0 and -old>max_gain:
                    max_gain=-old;witness={'subcell':i,'angle_cell':p,'row':name,'old_lower':str(old),'new_lower':str(selected),'new_pair':kbest,'closed_angle':[str(l),str(r)]}
                if selected<worstpose:worstpose,worstposeid=selected,[i,p,name,kbest]
    return {'reference_bin':j,'common_bracket':[str(left),str(right)],'parameter_cells':end-begin,'angle_cells_per_parameter_cell':angular,'root_margin_lower':str(worstroot),'root_id':worstrootid,'placement_margin_lower':str(worstpose),'placement_id':worstposeid,'compensated_comparisons':changed,'example':witness}


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());extra=[tuple(map(F,r)) for r in cfg['lines']];ns,src=geometry(archive,extra);results=[]
    for jlo,jhi,l,r in cfg['classes']:
        for j in range(jlo,jhi):
            result=root_conditioned(ns,j,F(l),F(r),cfg['angle_cells']);results.append(result);print(j,result['root_margin_lower'],result['placement_margin_lower'],result['compensated_comparisons'],flush=True)
    assert all(arb(v['root_margin_lower'])>0 and arb(v['placement_margin_lower'])>0 for v in results)
    out.write_text(json.dumps({'input':cfg,'results':results,'geometry_sha256':sha256(src.encode()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'same_anchor_kernel_sha256':sha256(Path(__file__).with_name('certify_ah_coupled_brackets.py').read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'bits':128},indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
