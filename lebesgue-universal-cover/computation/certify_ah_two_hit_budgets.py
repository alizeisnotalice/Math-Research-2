"""Couple first and second cap-depth budgets, and detect when the second hit must be in the old ear."""
from fractions import Fraction as F
import contextlib,inspect,io
from certify_ah_first_hit_chain import geometry as first_geometry
from certify_ah_single_contact_chain import support_bound as base_support
from certify_ah_attainable_root import validate as base_validate
from certify_ah_margin_envelope import replace_line


def geometry(archive,extra,*,far,negative,eps):
    old,prefix,_=first_geometry(archive,extra,far=far,negative=negative,eps=eps)
    prefix=replace_line(prefix,'parameter_boxes=[]','physical_deficits=[]\nparameter_boxes=[]')
    prefix=replace_line(prefix,'    parameter_boxes.append((l,pl,used,a))','    physical_deficits.append((l,cell_deficit))\n    parameter_boxes.append((l,pl,used,a))')
    ns={'__file__':old['__file__'],'__name__':'__main__'}
    with contextlib.redirect_stdout(io.StringIO()):exec(compile(prefix,old['__file__'],'exec'),ns)
    ns['two_hit_budgets']={}
    for j,E in enumerate(ns['deficit_by_bin']):
        if E.is_zero():continue
        A=ns['A'];a,b=ns['knots'][j:j+2];l1,l2=[ns['chart'](ns['add'](ns['z'],ns['mul'](A(w),ns['D'])))[0] for w in (a,b)];l=l1.union(l2);lp=ns['coord_hull'](ns['psi'],l);tab=[]
        for q,center in enumerate((lp,l)):
            first=ns['labelled_deficits'][q][j];partner=center+(-A(F(2,3))*first).union(A(0));second=max([A(0)]+[ee for ell,ee in ns['physical_deficits'] if ell.overlaps(partner)]);tab.append((first,second))
        ns['two_hit_budgets'][j]=tab
    return ns,prefix,None


def support_bound(ns,j,a,b):
    template=inspect.getsource(base_support);mark='    # Only the coupled bound is exported.'
    fragment='''    src=replace_line(src,'def value(x,eta,xi,first,de,partner,m):','def value(x,eta,xi,first,de,partner,m,orientation):')
    former="return first+de-partner+support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de,-first),zero,A(2))"
    new="""kappa=support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de),zero,A(1))
            first_bound,second_bound=ns['two_hit_budgets'][j][orientation]
            e0=min(first_bound,second_bound)
            def at(e,v):return ((kappa-first)*e+(de-partner+kappa)*v)/E
            return extremal((zero,at(e0,e0),at(first_bound,zero),at(first_bound,e0)))"""
    src=replace_line(src,former,new)
    src=replace_line(src,'value(ll,lr,a2,t4,dt,t2,local_m[0]),value(lr,ll,a4,t2,-dt,t4,local_m[1])','value(ll,lr,a2,t4,dt,t2,local_m[0],0),value(lr,ll,a4,t2,-dt,t4,local_m[1],1)')
'''
    template=replace_line(template,mark,fragment+mark);sc=dict(base_support.__globals__);exec(compile(template,'two_hit_budget_support','exec'),sc);return sc['support_bound'](ns,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(base_validate);sc={**base_validate.__globals__,'support_bound':support_bound};exec(compile(src,'two_hit_filtered_root','exec'),sc);return sc['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
