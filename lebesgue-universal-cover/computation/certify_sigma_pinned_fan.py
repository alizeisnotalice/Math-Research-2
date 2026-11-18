"""Subdivide the ideal vertex fans by the six exactly pinned hexagon supports."""
import inspect
from certify_ah_same_fan_cancellation import corrected_support
from certify_ah_two_hit_budgets import support_bound as budget_bound
from certify_ah_local_physical_support import support_bound as local_support
from certify_ah_margin_envelope import replace_line


def build(ns,edges):
    cross=ns['cross'];norms=list(map(ns['av'],ns['d']['n']));fans=[]
    for k,(vertex,arc) in enumerate(ns['sector']):
        if arc:continue
        t=(k+1)%10;e,f=edges[k],edges[t];pins=[u for u in norms if cross(e,u)>0 and cross(u,f)>0]
        pins.sort(key=lambda u:float(cross(e,u)/ns['dot'](e,u)) if ns['dot'](e,u)>0 else 100.0)
        rays=[(e,ns['orders'][k][1])]+[(u,None) for u in pins]+[(f,ns['orders'][t][1])]
        for (l,i),(r,j) in zip(rays,rays[1:]):
            assert cross(l,r)>0
            fans.append((vertex,l,r,i,j))
    return fans


def joint(ns,j,a,b):
    src=inspect.getsource(budget_bound)
    mark='    # Only the coupled bound is exported.'
    start="    def choices(c,u):";stop="    def extremal(vals):"
    part='''    fans=ns['_build_pinned_fans'](ns,edges)
    def pieces(c,u):
        groups=[]
        for sign in (1,-1):
            g=mul(sign,u);coeff=sign*c
            for idx,(vertex,l,r,il,ir) in enumerate(fans):
                if cross(l,g)<0 or cross(g,r)<0:continue
                den=cross(l,r);ll=cross(g,r)/den;rr=cross(l,g)/den
                groups.append((idx,coeff,g,ll,rr))
        assert groups
        return groups
    def choices(c,u):
        actions=[]
        for idx,coeff,g,ll,rr in pieces(c,u):
            vertex,l,r,il,ir=fans[idx]
            if not coeff>0:actions.append([(vertex,mul(coeff,g))])
            if not coeff<0:actions.append([(i,mul(v,ray)) for i,v,ray in ((il,coeff*ll,l),(ir,coeff*rr,r)) if i is not None])
        return actions
    def correction(terms):
        ledger={};cor=[]
        for c,u in terms:
            modes=pieces(c,u)
            if len(modes)!=1:continue
            idx,w,g,ll,rr=modes[0]
            if not ll>0 or not rr>0:continue
            pos,neg=ledger.setdefault(idx,([],[]))
            if w>0:pos.append([w*ll,w*rr])
            if w<0:neg.append([-w,ll,rr,g])
        for idx,(pos,neg) in ledger.items():
            vertex,l,r,il,ir=fans[idx]
            for demand,ll,rr,g in neg:
                for avail in pos:
                    take=min(demand.lower(),(avail[0].lower()/ll.upper()).lower(),(avail[1].lower()/rr.upper()).lower())
                    if not take>0:continue
                    if il is not None:cor.append((il,mul(-take*ll,l)))
                    if ir is not None:cor.append((ir,mul(-take*rr,r)))
                    cor.append((vertex,mul(take,g)));avail[0]-=take*ll;avail[1]-=take*rr;demand-=take
                    if not demand>0:break
        return cor
'''
    command=f"    src=src[:src.index({start!r})]+{part!r}+src[src.index({stop!r}):]\n"
    old='common=extremal([event(sum(g,[])) for g in itertools.product(*groups)])'
    new='adjust=correction(terms);common=extremal([event(sum(g,[])+adjust) for g in itertools.product(*groups)])'
    command+='    src=replace_line(src,'+repr(old)+','+repr(new)+')\n'
    anchor='template=replace_line(template,mark,fragment+mark)';src=replace_line(src,anchor,'fragment+=command;'+anchor)
    if ns.get('_sharpen_jets'):
        extra="""    narrow={}
    for ori in range(2):
        Fh,Sh=ns['two_hit_budgets'][j][ori];vh=min(Fh,Sh)
        for k,ell,T in ((4,l,t),(2,lp,-t)):
            is_first=(k==4)==(ori==0)
            dl=zero if is_first else (-A(F(2,3))*Fh).union(zero)
            dt=(-Fh).union(zero) if is_first else (-vh).union(Fh)
            su=z;sv=z;parts=8
            for rj in range(parts):
                lam=A(F(rj,parts)).union(A(F(rj+1,parts)))
                x0,ll,nn=ns['chart_bounds'](ell+lam*dl,T+lam*dt);dx,dy=ns['derivs'](x0,ll,nn)
                su=add(su,mul(A(F(1,parts)),dx));sv=add(sv,mul(A(F(1,parts)),dy))
            narrow[ori,k]=(su,sv) if k==4 else (reflection(su),reflection(sv))
"""
        step='    src=replace_line(src,'+repr('    fans=ns[\'_build_pinned_fans\'](ns,edges)')+','+repr(extra+"    fans=ns['_build_pinned_fans'](ns,edges)")+')\n'
        orig='value(ll,lr,a2,t4,dt,t2,local_m[0],0),value(lr,ll,a4,t2,-dt,t4,local_m[1],1)'
        conv='value(ll,lr,dot(vec[2],narrow[0,2][0]),dot(vec[4],narrow[0,4][1]),dt,dot(vec[2],narrow[0,2][1]),local_m[0],0),value(lr,ll,dot(vec[4],narrow[1,4][0]),dot(vec[2],narrow[1,2][1]),-dt,dot(vec[4],narrow[1,4][1]),local_m[1],1)'
        command+=step+'    src=replace_line(src,'+repr(orig)+','+repr(conv)+')\n'
    if '_contact_levels' in ns:
        before='kappa=support_chain((x/(m*A(\'.95\')),-eta/m,-xi/m,partner,-de),zero,A(1))'
        after="""coeff=(x/(m*A('.95')),-eta/m,-xi/m,partner,-de)
            levels=ns['_contact_levels'];offset=0;kappa=zero
            for level,count in enumerate(levels):
                if count:kappa+=support_chain(coeff[offset:offset+count],A(F(level,len(levels))),A(F(level+1,len(levels))))
                offset+=count
            assert offset==5"""
        src=replace_line(src,before,after)
        src=replace_line(src,'at(e0,e0),at(first_bound,zero),at(first_bound,e0)','at(e0,e0),at(first_bound,e0)')
    ns['_build_pinned_fans']=build
    gs={**budget_bound.__globals__,'command':command}
    exec(compile(src,'pinned_fan_budget','exec'),gs);return gs['support_bound'](ns,j,a,b)


def support_bound(ns,j,a,b):
    gs={**local_support.__globals__,'base_support':joint};exec(compile(inspect.getsource(local_support),'pinned_local_budget','exec'),gs);return gs['support_bound'](ns,j,a,b)


def combined(ns,j,a,b):
    from certify_ah_same_fan_cancellation import support_bound as old
    E,v,edges,fn,_=support_bound(ns,j,a,b);other=old(ns,j,a,b)[3]
    def bound(actions):return fn(actions).min(other(actions))
    return E,v,edges,bound,None


def conditional(*args,**kw):
    from certify_ah_affine_root_class import validate
    env={**validate.__globals__,'support_bound':combined}
    exec(compile(inspect.getsource(validate),'pinned_fan_root_choice','exec'),env);return env['validate'](*args,**kw)
