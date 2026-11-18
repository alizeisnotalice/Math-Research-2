"""Rigorous obstruction from a deformed Reuleaux pentagon with all non-dihedral compatibility roots excluded."""
import argparse,json,sys
from hashlib import sha256
from pathlib import Path
from fractions import Fraction as F
from flint import arb
from certify_sigma_extended_graph import geometry


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());ns,_,_=geometry(archive,[(F('.429'),F('.009'))],sigma=cfg['sigma'],far=F('.273'),negative=F('.306'),eps=F('.0102'));A,add,sub,mul,dot,norm,perp,qa=(ns[k] for k in ('A','add','sub','mul','dot','norm','perp','reflection'));half=A('.5');zero=A(0);w,e,p,d=map(F,[cfg[k] for k in ('w','deformation','pivot','angle')]);r=ns['r3'];normals=list(map(ns['av'],ns['d']['n']));mp=list(map(ns['av'],ns['d']['m']));ell,t=ns['chart'](add(ns['z'],mul(A(w),ns['D'])));partner=ns['psi'](ell);t-=A(e)
    L=add(ns['l0'],mul(ell,ns['el']));M=add(ns['l0'],mul(partner,ns['el']));N=add(ns['mid'],mul(t,ns['en']));V=[L,qa(M),qa(ns['X'](partner,-t)),N,ns['X'](ell,t)];edges=[sub(V[i],V[j]) for i,j in ns['orders']];orders=ns['orders'];cross=ns['cross'];pent=min(A(1)-dot(sub(V[k],V[l]),sub(V[k],V[l])) for k in range(5) for l in range(k) if (k,l) not in orders);assert pent>0
    assert all(dot(vec,vec).contains(1) for vec in edges) and all(cross(x,edges[(i+1)%10])>0 for i,x in enumerate(edges));trig=ns['trig'];sup=ns['supports'];F_at=ns['F_at'];n0=normals[0]
    pin=[]
    for k,i in enumerate((1,3,0)):
        u=mul(-1,normals[k]);m=[j for j,vv in enumerate(edges) if cross(vv,u)>0 and cross(u,edges[(j+1)%10])>0];assert len(m)==1 and tuple(ns['sector'][m[0]])==(i,0),(k,m);assert dot(V[i],u).contains(half);pin.append(m[0])
    S=sum((-1)**k*dot(V[i],perp(normals[k])) for k,i in enumerate((1,3,0)));sg=1 if S>0 else -1;endpoint=F('.10');first=(sg*S-A('1.5')*(A(endpoint/2)*ns['angle_unit']).tan()).lower();assert first>0,(S,first)
    low=A(1);cells=0;deepest=0
    def visit(l,r,depth):
        nonlocal low,cells,deepest
        m=(l+r)/2;fm=F_at(V,edges,m);dw=ns['box'](l-m,r-m)*ns['angle_unit'];fp=sum((-1)**k*sup(V,edges,trig(ns['box'](l+60*k,r+60*k)),n0,True) for k in range(3));g=(sg*(fm+dw*fp)).lower()
        if g>0:low=min(low,g);cells+=1;deepest=max(deepest,depth);return
        assert depth<20,(l,r,fm,g,sg);visit(l,m,depth+1);visit(m,r,depth+1)
    for k in range(240):visit(endpoint+(60-2*endpoint)*k/240,endpoint+(60-2*endpoint)*(k+1)/240,0)
    direct=add(mul(A(d).cos(),ns['direct']),mul(A(d).sin(),ns['normal']));cu=mul(-1,perp(direct));pivot=add(ns['z'],mul(A(p),ns['D']));rows=[('m3',mp[3],half),('m5',mp[5],half),('AH-cut',cu,dot(cu,pivot)),('AH-original',mul(-1,ns['normal']),-dot(ns['normal'],ns['z']))]
    if cfg.get('cut_only'):rows=rows[:3]
    for k in range(0 if cfg.get('cut_only') else 161):
        z=F(k,160);q=ns['q'](A(z));dq=qa(ns['point_jets'](A(z))[1][2]);u=mul(-1/norm(dq),perp(dq));rows.append((f'paired-tangent:{z}',u,dot(u,q)))
    bounds=[]
    for turn in range(0,360,60):
        for sign in [1,-1]:
            val,label=max(((sup(V,edges,trig(A(-turn)),(u[0],sign*u[1]))-c).lower(),k) for k,u,c in rows);assert val>0,(turn,sign,label,val);bounds.append({'turn':turn,'sign':sign,'edge':label,'violation':str(val)})
    obj={'parameters':cfg,'five_vertices':[[str(z) for z in point] for point in V],'nonedge_diameter_squared_gap_lower':str(pent.lower()),'pinned_vertex_sectors':pin,'sign_inside_0_to_60':sg,'endpoint_neighborhood_degrees':str(endpoint),'endpoint_signed_normalized_margin':str(first),'interior_signed_lower':str(low),'interior_cells':cells,'max_depth':deepest,'dihedral_witnesses':bounds,'all_pose_violation_lower':str(min(arb(s['violation']).lower() for s in bounds)),'python':sys.version,'flint':ns['flint_version'],'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()};out.write_text(json.dumps(obj,indent=2)+'\n');print(cells,sg,first,low,obj['all_pose_violation_lower'],bounds)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('input',type=Path);p.add_argument('output',type=Path);o=p.parse_args();run(o.archive.resolve(),o.input,o.output)
