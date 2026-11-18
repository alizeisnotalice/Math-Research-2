"""Cancel opposite support coefficients within a common ideal vertex cone by true support subadditivity."""
import inspect
from certify_ah_margin_envelope import replace_line
from certify_ah_two_hit_budgets import support_bound as two_support
from certify_ah_local_physical_support import support_bound as physical_support
from certify_ah_affine_root_class import validate as affine_validate


CORRECTION='''    def shared_fan_correction(terms):
        by_vertex={}
        for c,u in terms:
            cr=[cross(e,u) for e in edges]
            active=[k for k in range(10) if not cr[k]<0 and not cr[(k+1)%10]>0]
            if len(active)!=1:continue
            k=active[0]
            if ns['sector'][k][1]:
                k=(k+5)%10;u=mul(-1,u);c=-c
            vert,arc=ns['sector'][k];assert not arc
            r=(k+1)%10;e0,e1=edges[k],edges[r];det=cross(e0,e1)
            ll,rr=cross(u,e1)/det,cross(e0,u)/det
            if ll<=0 or rr<=0 or not (c>0 or c<0):continue
            pos,neg=by_vertex.setdefault(k,([],[]))
            if c>0:pos.append([c*ll,c*rr])
            else:neg.append([-c,ll,rr,u])
        corrections=[]
        for k,(pos,neg) in by_vertex.items():
            rr=(k+1)%10;vert=ns['sector'][k][0];i0,i1=ns['orders'][k][1],ns['orders'][rr][1]
            for entry in neg:
                demand,ll,r,u=entry
                if not ll>0 or not r>0:continue
                for avail in pos:
                    take=min(demand.lower(),(avail[0].lower()/ll.upper()).lower(),(avail[1].lower()/r.upper()).lower())
                    if not take>0:continue
                    corrections.extend(((i0,mul(-take*ll,edges[k])),(i1,mul(-take*r,edges[rr])),(vert,mul(take,u))))
                    avail[0]-=take*ll;avail[1]-=take*r;demand-=take
                    if not demand>0:break
        return corrections
'''


def corrected_support(ns,j,a,b):
    src=inspect.getsource(two_support)
    command="    src=replace_line(src,"+repr('    def bound(terms,detail=False):')+","+repr(CORRECTION+'    def bound(terms,detail=False):')+")\n"
    former='common=extremal([event(sum(g,[])) for g in itertools.product(*groups)])'
    new='adjust=shared_fan_correction(terms); common=extremal([event(sum(g,[])+adjust) for g in itertools.product(*groups)])'
    # Both scalar bounds apply; retain the smaller when interval enclosures lose correlation.
    new+=';common=extremal([event(sum(g,[])) for g in itertools.product(*groups)]).min(common) if adjust else common'
    command+="    src=replace_line(src,"+repr(former)+","+repr(new)+")\n"
    marker='template=replace_line(template,mark,fragment+mark)'
    src=replace_line(src,marker,'fragment+=extra_patch;'+marker)
    env={**two_support.__globals__,'extra_patch':command}
    exec(compile(src,'subadditive_fan_events','exec'),env);return env['support_bound'](ns,j,a,b)


def support_bound(ns,j,a,b):
    env={**physical_support.__globals__,'base_support':corrected_support}
    exec(compile(inspect.getsource(physical_support),'common_fan_physical_subcells','exec'),env)
    return env['support_bound'](ns,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=4,anchor=None,threshold=None,condition_sign=1):
    env={**affine_validate.__globals__,'support_bound':support_bound}
    exec(compile(inspect.getsource(affine_validate),'shared_fan_signed_pose','exec'),env)
    return env['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth,anchor,threshold,condition_sign)
