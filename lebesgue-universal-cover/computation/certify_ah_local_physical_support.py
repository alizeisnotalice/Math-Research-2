"""Intersect each fine reference cell with the true first hit and its possible partner in physical chart coordinates."""
from bisect import bisect_left,bisect_right
from fractions import Fraction as F
import inspect
from certify_ah_two_hit_budgets import support_bound as base_support
from certify_ah_pinned_normals import validate as pinned_validate


def local_budgets(ns,a,b):
    A=ns['A']
    if '_physical_lookup' not in ns:
        rows=sorted(ns['physical_deficits'],key=lambda row:row[0].lower())
        lowers=[ell.lower() for ell,_ in rows]
        prefix=[];top=rows[0][0].upper()
        for ell,_ in rows:
            top=max(top,ell.upper());prefix.append(top)
        ns['_physical_lookup']=rows,lowers,prefix
    rows,lowers,prefix=ns['_physical_lookup']
    def lookup(target):
        start=bisect_left(prefix,target.lower());end=bisect_right(lowers,target.upper())
        return max([A(0)]+[d for ell,d in rows[start:end] if ell.overlaps(target)])
    coords=[ns['chart'](ns['add'](ns['z'],ns['mul'](A(x),ns['D'])))[0] for x in (a,b)]
    l=coords[0].union(coords[1]);p=ns['coord_hull'](ns['psi'],l)
    budgets=[]
    for first,other in ((l,p),(p,l)):
        e=lookup(first);v=lookup(other+(-A(F(2,3))*e).union(A(0)))
        budgets.append((e,v))
    return budgets


def support_bound(ns,j,a,b):
    new=local_budgets(ns,a,b);old=ns['two_hit_budgets'][j]
    budgets=[(min(e0,e),min(v0,v)) for (e0,v0),(e,v) in zip(old,new,strict=True)]
    if '_first_orientation' in ns:
        budgets[1-ns['_first_orientation']]=(ns['A'](0),ns['A'](0))
    # The old normalization (including its nonzero denominator) remains valid.
    local=dict(ns);local['two_hit_budgets']={j:budgets}
    return base_support(local,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(pinned_validate)
    old_line='a=F(j,100)+F(i,7500);b=a+F(1,7500);'
    new_line="ref_lo,ref_hi=ns['knots'][j:j+2];step=(ref_hi-ref_lo)/75;a=ref_lo+i*step;b=a+step;"
    mark='src=inspect.getsource(base_validate)'
    src=src.replace(mark,mark+';src=replace_line(src,'+repr(old_line)+','+repr(new_line)+')',1)
    src=src.replace("**pose_validate.__globals__,'pinned_source':src","**pose_validate.__globals__,'support_bound':support_bound,'pinned_source':src")
    env={**pinned_validate.__globals__,'support_bound':support_bound}
    exec(compile(src,'local_physical_support_rows','exec'),env)
    return env['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
