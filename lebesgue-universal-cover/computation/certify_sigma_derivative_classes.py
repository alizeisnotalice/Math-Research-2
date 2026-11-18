"""Condition root and placement bounds on the compatibility derivative of the same body at zero."""
from fractions import Fraction as F
import inspect,itertools
from certify_sigma_pinned_fan import support_bound
from certify_ah_affine_root_class import validate as affine
from certify_ah_margin_envelope import replace_line
from certify_ah_pinned_normals import pinned_equation


def derivative_bound(ns,j,a,b):
    A=ns['A'];E,v,ed,base,_=support_bound(ns,j,a,b);env=inspect.getclosurevars(base).nonlocals;choices,correction,extremal,event=[env[k] for k in ('choices','correction','extremal','event')];normals=list(map(ns['av'],ns['d']['n']));actions=[(i,ns['mul'](A((-1)**k),(-normals[k][1],normals[k][0]))) for k,i in enumerate([1,3,0])];S0=sum(ns['dot'](v[i],u) for i,u in actions);ns['_model_slope']=S0
    def bound(ts):
        extra=[];terms=[]
        for c,u in ts:
            if u is None:extra.extend((i,ns['mul'](c,w)) for i,w in actions)
            else:terms.append((c,u))
        if not extra:return base(terms)
        cor=correction(terms);groups=[choices(c,u) for c,u in terms];return extremal([event(sum(g,[])+cor+extra) for g in itertools.product(*groups)]).upper()
    ns['_actual_slope']=(S0-E*bound([(A(-1),None)])).lower().union((S0+E*bound([(A(1),None)])).upper());return E,v,ed,bound,None


def validate(ns,j,left,right,turn,sign,orientation,*,cut,side,begin=0,end=75,angular=12,depth=4):
    src=inspect.getsource(affine)
    src=replace_line(src,'assert orientation==1 and condition_sign in (-1,1)','assert orientation in (1,-1) and condition_sign in (-1,1)')
    src=replace_line(src,'anchor_terms=[(condition_sign*c,u) for c,u in root_terms(A(anchor))]','anchor_terms=[(A(condition_sign),None)]')
    src=replace_line(src,"f_anchor=condition_sign*(ns['F_at'](v,e,anchor)+A(threshold))","f_anchor=condition_sign*(ns['_model_slope']-A(threshold))")
    src=replace_line(src,'for phi,s in ((left,1),(right,-1)):','for phi,s in ((left,orientation),(right,-orientation)):')
    src=replace_line(src,'if phi==anchor and s==condition_sign:slack=max(slack,-condition_sign*A(threshold))',"S=ns['_actual_slope']; clipped=(S.max(A(threshold)).lower().union(S.upper()) if condition_sign==1 else S.lower().union(S.min(A(threshold)).upper()));slack=max(slack,(s*pinned_equation(ns,clipped,A(phi))).lower())")
    src=replace_line(src,'if gap>0:',"theta=interval*ns['angle_unit'];rad=A(F(3,2))*(theta/2).tan();pin=theta.sin()*(clipped+(-rad).union(rad));gap=max(gap,pin.lower(),(-pin).lower())\n            if gap>0:")
    src=replace_line(src,'F(1),F(11,10)', 'F(1,100),F(1,50),F(1,25),F(1,16),F(3,32),F(1,8),F(3,16),F(1,4),F(1),F(11,10)')
    env={**affine.__globals__,'support_bound':derivative_bound,'pinned_equation':pinned_equation};exec(compile(src,'split_on_actual_initial_slope','exec'),env);return env['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth,anchor=F(0),threshold=cut,condition_sign=side)
