"""A second certified root under the same body's failure of the first sign test."""
from fractions import Fraction as F
from certify_ah_local_physical_support import support_bound


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=4,anchor=F(23),threshold=F('0.00003'),condition_sign=1):
    assert orientation==1 and condition_sign in (-1,1)
    A=ns['A'];trig,rot,sup,cross=[ns[k] for k in ('trig','rotate','supports','cross')]
    zero=A(0);one=A(1);half=A(F(1,2));normal=list(map(ns['av'],ns['d']['n']));n0=normal[0]
    def root_terms(angle):
        return [(A((-1)**k),rot(n0,trig(angle+A(60*k)))) for k in range(3)]
    anchor_terms=[(condition_sign*c,u) for c,u in root_terms(A(anchor))];host=[]
    for label,u,c in ns['rows']:
        options=[]
        for k in range(3):
            n,m=normal[k:k+2];det=cross(n,m);options.append((k,cross(u,m)/det,cross(n,u)/det))
        host.append((label,u,c,options))
    poses=[(turn,sign)]+[(k,s) for k in (0,60,120,180,240,300) for s in (1,-1) if (k,s)!=(turn,sign)]
    rm=pm=one;rid=pid=None;excluded=placed=empty=0;counts={};gap_min=one
    for i in range(begin,end):
        r0,r1=ns['knots'][j:j+2];step=(r1-r0)/75;a=r0+i*step;b=a+step
        E,v,e,bound,_=support_bound(ns,j,a,b)
        f_anchor=condition_sign*(ns['F_at'](v,e,anchor)+A(threshold))
        if f_anchor+E*bound(anchor_terms)<0:
            empty+=1;continue
        def lower(model,minus_terms,target=zero):
            value=(model-E*bound(minus_terms)).lower()
            if value>target:return value
            for q in (F(1),F(11,10),F(21,20),F(23,20),F(5,4),F(7,8),F(1,2),F(2),F(1,4),F(3,2),F(3,4),F(3),F(1,8)):
                g=A(q);mix=minus_terms+[(g*c,u) for c,u in anchor_terms]
                value=max(value,(model-g*f_anchor-E*bound(mix)).lower())
                if value>target:break
            return value
        for phi,s in ((left,1),(right,-1)):
            terms=root_terms(A(phi));slack=lower(s*ns['F_at'](v,e,phi),[(-s*c,u) for c,u in terms])
            if phi==anchor and s==condition_sign:slack=max(slack,-condition_sign*A(threshold))
            if slack<rm:rm,rid=slack,[i,str(phi)]
        def visit(l,r,lev):
            nonlocal pm,pid,placed,excluded,gap_min
            mid=(l+r)/2;interval=ns['box'](l,r);dw=ns['angle_unit']*ns['box'](l-mid,r-mid);fm=ns['F_at'](v,e,mid)
            deriv=sum((-1)**k*sup(v,e,trig(interval+A(60*k)),n0,True) for k in range(3));fm+=dw*deriv
            terms=root_terms(interval);gap=lower(fm,[(-c,u) for c,u in terms])
            if gap<=0:gap=lower(-fm,terms)
            if gap>0:
                excluded+=1;gap_min=min(gap_min,gap);return
            best=None
            for t,s in poses:
                mt,sn=trig(A(mid-t)),trig(interval-A(t));pre=[(u[0],s*u[1]) for u in normal]
                h=[sup(v,e,mt,n) for n in pre];d=[sup(v,e,sn,n,True) for n in pre];U=[rot(n,sn) for n in pre]
                worst=one;edge=None
                for label,u,c,opts in host:
                    x=(u[0],s*u[1]);X=rot(x,sn);hh=sup(v,e,mt,x);dh=sup(v,e,sn,x,True)
                    margin=None
                    for k,al,be in opts:
                        model=c-hh+al*(h[k]-half)+be*(h[k+1]-half)+dw*(-dh+al*d[k]+be*d[k+1])
                        value=(model-A('1.07')*(one+abs(al)+abs(be))*E).lower()
                        if value<A('.0001'):value=lower(model,[(one,X),(-al,U[k]),(-be,U[k+1])])
                        if margin is None or value>margin:margin,kbest=value,k
                        if margin>A('.0001'):break
                    if margin<=0:
                        _,al0,be0=opts[0]
                        for mu in (F(1,4),-F(1,4),F(1,2),-F(1,2),F(3,4),-F(3,4),F(1),-F(1),F(1,8),-F(1,8),F(3,2),-F(3,2),F(2),-F(2)):
                            coeffs=(al0+A(mu),be0-A(mu),A(mu))
                            model=c-hh+sum(cc*(hz-half) for cc,hz in zip(coeffs,h[:3]))+dw*(-dh+sum(cc*dz for cc,dz in zip(coeffs,d[:3])))
                            value=lower(model,[(one,X)]+[(-cc,uz) for cc,uz in zip(coeffs,U[:3])])
                            if value>margin:margin,kbest=value,str(mu)
                            if margin>0:break
                    if margin<worst:worst,edge=margin,(label,kbest)
                    if margin<=0:break
                if best is None or worst>best[0]:best=(worst,edge,t,s)
                if worst>0:break
            val,edge,t,s=best
            if val<=0 and lev<depth:
                visit(l,mid,lev+1);visit(mid,r,lev+1);return
            placed+=1;key=f'{t},{s}';counts[key]=counts.get(key,0)+1
            if val<pm:pm,pid=val,[i,str(l),str(r),edge,t,s]
        for k in range(angular):visit(left+(right-left)*F(k,angular),left+(right-left)*F(k+1,angular),0)
    return {'reference_bin':j,'subcells':[begin,end],'bracket':[str(left),str(right)],'anchor_degrees':str(anchor),'anchor_comparison_sign':condition_sign,'anchor_threshold':str(-threshold),'root_margin_lower':str(rm),'root_id':rid,'placement_margin_lower':str(pm),'placement_id':pid,'root_free_leaves':excluded,'sign_separation_lower':str(gap_min),'placed_leaves':placed,'pose_counts':counts,'empty_conditional_subcells':empty}
