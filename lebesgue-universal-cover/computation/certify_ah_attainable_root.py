"""Reject uniformly root-free angles before testing placement on the remaining closed angular cells."""
from fractions import Fraction as F
from certify_ah_coupled_brackets import support_bound


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    A=ns['A'];cross=ns['cross'];rot=ns['rotate'];trig=ns['trig'];sup=ns['supports'];one=A(1);half=A('.5');zero=A(0)
    normals=list(map(ns['av'],ns['d']['n']));pre=[(v[0],sign*v[1]) for v in normals];n0=normals[0]
    rows=[]
    for name,u,c in ns['rows']:
        options=[]
        for k in range(3):
            n,m=normals[k:k+2];det=cross(n,m);assert det>0;alpha,beta=cross(u,m)/det,cross(n,u)/det;options.append((k,alpha,beta))
        rows.append((name,(u[0],sign*u[1]),c,options))
    rootmin=posemin=one;root_id=pose_id=None;leaves=reject=0;separation=one
    for i in range(begin,end):
        a=F(j,100)+F(i,7500);b=a+F(1,7500);E,v,e,bound,_=support_bound(ns,j,a,b)
        def root_terms(angle):return [(A((-1)**k),rot(n0,trig(angle+A(60*k)))) for k in range(3)]
        for phi,s in ((left,orientation),(right,-orientation)):
            terms=root_terms(A(phi));slack=(s*ns['F_at'](v,e,phi)-E*bound([(-s*c,u) for c,u in terms])).lower()
            if slack<rootmin:rootmin,root_id=slack,[i,str(phi)]
        def inspect(l,r,lev):
            nonlocal leaves,reject,posemin,pose_id,separation
            mid=(l+r)/2;interval=ns['box'](l,r); dw=ns['angle_unit']*ns['box']((l-r)/2,(r-l)/2);fmid=ns['F_at'](v,e,mid)
            derivative=sum((-1)**k*sup(v,e,trig(interval+A(60*k)),n0,True) for k in range(3));model=fmid+dw*derivative;terms=root_terms(interval);gap=model-E*bound([(-c,u) for c,u in terms]);
            no_root=gap.lower()
            if no_root<=0:no_root=(-model-E*bound(terms)).lower()
            if no_root>0:
                reject+=1;separation=min(separation,no_root);return
            mt=trig(A(mid-turn));sn=trig(interval-A(turn));hm=[sup(v,e,mt,n) for n in pre];dh=[sup(v,e,sn,n,True) for n in pre];un=[rot(n,sn) for n in pre];local=[]
            for name,u,c,options in rows:
                um=sup(v,e,mt,u);du=sup(v,e,sn,u,True);ru=rot(u,sn);low=None
                for k,alpha,beta in options:
                    i1=k+1;model=c-um+alpha*(hm[k]-half)+beta*(hm[i1]-half)+dw*(-du+alpha*dh[k]+beta*dh[i1]);rough=(model-A('1.07')*(one+abs(alpha)+abs(beta))*E).lower()
                    val=rough if rough>A('.0001') else (model-E*bound([(one,ru),(-alpha,un[k]),(-beta,un[i1])])).lower()
                    if low is None or val>low:low,kbest=val,k
                    if low>A('.0001'):break
                if low<=0 and lev<depth:
                    inspect(l,mid,lev+1);inspect(mid,r,lev+1);return
                local.append((low,name,kbest))
            leaves+=1
            low,name,kbest=min(local,key=lambda r:r[0]);
            if low<posemin:posemin,pose_id=low,[i,str(l),str(r),name,kbest]
        for k in range(angular):inspect(left+(right-left)*F(k,angular),left+(right-left)*F(k+1,angular),0)
    return {'reference_bin':j,'subcells':[begin,end],'bracket':[str(left),str(right)],'pose':[turn,sign,orientation],'root_margin_lower':str(rootmin),'root_id':root_id,'placement_margin_lower':str(posemin),'placement_id':pose_id,'placed_leaves':leaves,'root_free_leaves':reject,'sign_separation_lower':str(separation) if reject else None,'starting_angle_cells':angular,'max_subdivision_depth':depth}
