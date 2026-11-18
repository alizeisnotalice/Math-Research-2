"""At one genuine compatibility root, verify one common dihedral placement per angular leaf."""
import inspect
from certify_ah_attainable_root import validate as base_validate
from certify_ah_two_hit_budgets import support_bound
from certify_ah_margin_envelope import replace_line


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(base_validate)
    src=replace_line(src,"rows.append((name,(u[0],sign*u[1]),c,options))","rows.append((name,u,c,options))")
    src=replace_line(src,'rootmin=posemin=one;', 'pose_count={}\n    rootmin=posemin=one;')
    start=src.index('            mt=trig(A(mid-turn));');endpos=src.index('        for k in range(angular):',start)
    new='''            best=None
            available=[(turn,sign)]+[(t,s) for t,s in ((300,1),(120,-1),(0,1),(0,-1),(60,1),(60,-1),(120,1),(180,1),(180,-1),(240,1),(240,-1),(300,-1)) if (t,s)!=(turn,sign)]
            for t,s in available:
                mt=trig(A(mid-t));sn=trig(interval-A(t));pre=[(u[0],s*u[1]) for u in normals]
                hm=[sup(v,e,mt,n) for n in pre];dh=[sup(v,e,sn,n,True) for n in pre];un=[rot(n,sn) for n in pre];local=[]
                for name,u,c,options in rows:
                    uref=(u[0],s*u[1]);um=sup(v,e,mt,uref);du=sup(v,e,sn,uref,True);ru=rot(uref,sn);low=None
                    for k,alpha,beta in options:
                        i1=k+1;model=c-um+alpha*(hm[k]-half)+beta*(hm[i1]-half)+dw*(-du+alpha*dh[k]+beta*dh[i1]);rough=(model-A('1.07')*(one+abs(alpha)+abs(beta))*E).lower()
                        val=rough if rough>A('.0001') else (model-E*bound([(one,ru),(-alpha,un[k]),(-beta,un[i1])])).lower()
                        if low is None or val>low:low,kbest=val,k
                        if low>A('.0001'):break
                    local.append((low,name,kbest))
                    if low<=0:break
                worst=min(local,key=lambda r:r[0]);
                if best is None or worst[0]>best[0]:best=(*worst,t,s)
                if worst[0]>0:break
            low,name,kbest,t,s=best
            if low<=0 and lev<depth:
                inspect(l,mid,lev+1);inspect(mid,r,lev+1);return
            leaves+=1
            key=f'{t},{s}';pose_count[key]=pose_count.get(key,0)+1
            if low<posemin:posemin,pose_id=low,[i,str(l),str(r),name,kbest,t,s]
'''
    src=src[:start]+new+src[endpos:];src=replace_line(src,"'placed_leaves':leaves,","'placed_leaves':leaves,'verified_dihedral_pose_counts':pose_count,")
    sc={**base_validate.__globals__,'support_bound':support_bound};exec(compile(src,'root_dihedral_choices','exec'),sc);return sc['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
