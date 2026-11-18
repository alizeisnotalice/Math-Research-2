"""Couple the DC and DE parameter chains using the second hit's graph lower bound."""
import inspect
from fractions import Fraction as F
from certify_ah_coupled_brackets import support_bound as old_bound
from certify_ah_attainable_root import validate as old_validate
from certify_ah_margin_envelope import replace_line


def slopes(ns):
    if '_contact_chain_slopes' in ns:return ns['_contact_chain_slopes']
    A=ns['A'];rows=[]
    for k in range(430):
        a,b=F(k,1000),F(k+1,1000)
        ell=[ns['chart'](ns['add'](ns['z'],ns['mul'](A(x),ns['D'])))[0] for x in (a,b)];l=ell[0].union(ell[1]);p=ns['coord_hull'](ns['psi'],l)
        _,gl,gn=ns['points'](ns['box'](a,b));m=-ns['dot'](gn,ns['D'])/ns['dot'](gl,ns['D']);s=-ns['psi_prime'](ns['add'](ns['l0'],ns['mul'](l,ns['el'])));assert m>A('3.5') and A('.95')<s<A('1.05');rows.extend(((l,m.lower()),(p,(m/s).lower())))
    ns['_contact_chain_slopes']=rows;return rows


def local_slope(ns,interval,E):
    tube=interval+(-ns['A'](F(2,3))*E).union(ns['A'](0));values=[m for domain,m in slopes(ns) if domain.overlaps(tube)];assert values;return min(values)


def support_bound(ns,j,a,b):
    src=inspect.getsource(old_bound)
    src=replace_line(src,'    def choices(c,u):','    local_m=(local_slope(ns,lp,E),local_slope(ns,l,E))\n    def choices(c,u):')
    src=replace_line(src,'def value(x,eta,xi,first,de,partner):','def value(x,eta,xi,first,de,partner,m):')
    src=replace_line(src,"return support_chain((x/A('.95'),-eta,-xi),zero,A(F(2,3)))+support_chain((first,de,-partner),A(-1),A(1))", "return first+de-partner+support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de,-first),zero,A(2))")
    src=replace_line(src,'value(ll,lr,a2,t4,dt,t2),value(lr,ll,a4,t2,-dt,t4)','value(ll,lr,a2,t4,dt,t2,local_m[0]),value(lr,ll,a4,t2,-dt,t4,local_m[1])')
    # Only the coupled bound is exported. The two-simplex example enumerator belongs to its older model.
    src=src.split('    def incompatibility(')[0]+'    return E,v,edges,bound,None\n';sc={**old_bound.__globals__,'local_slope':local_slope};exec(compile(src,'single_contact_polytope','exec'),sc);return sc['support_bound'](ns,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(old_validate);sc={**old_validate.__globals__,'support_bound':support_bound};exec(compile(src,'feasible_angle_single_contact_chain','exec'),sc);return sc['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
