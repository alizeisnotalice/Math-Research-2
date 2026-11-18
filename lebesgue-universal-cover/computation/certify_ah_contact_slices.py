"""Partition the ordered common contact simplex; each class may choose its own genuine root bracket."""
from fractions import Fraction as F
import inspect
from certify_ah_two_hit_budgets import support_bound as two_support
from certify_ah_local_physical_support import support_bound as local_support,validate as local_validate
from certify_ah_margin_envelope import replace_line


def sliced_support(ns,j,a,b):
    src=inspect.getsource(two_support)
    original="kappa=support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de),zero,A(1))"
    changed="""coefficients=(x/(m*A('.95')),-eta/m,-xi/m,partner,-de)
            bins=ns.get('_contact_levels',(5,))
            kappa=zero;offset=0;denom=len(bins)
            for level,count in enumerate(bins):
                if count:kappa+=support_chain(coefficients[offset:offset+count],A(F(level,denom)),A(F(level+1,denom)))
                offset+=count
            assert offset==5"""
    src=replace_line(src,original,changed)
    # v is defined to be min(e,S) before any subdivision.  Its graph has two linear pieces.
    src=replace_line(src,'at(e0,e0),at(first_bound,zero),at(first_bound,e0)','at(e0,e0),at(first_bound,e0)')
    env=dict(two_support.__globals__);exec(compile(src,'sliced_contact_chain','exec'),env)
    return env['support_bound'](ns,j,a,b)


def support_bound(ns,j,a,b):
    env={**local_support.__globals__,'base_support':sliced_support}
    exec(compile(inspect.getsource(local_support),'sliced_physical_cells','exec'),env)
    return env['support_bound'](ns,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    env={**local_validate.__globals__,'support_bound':support_bound}
    exec(compile(inspect.getsource(local_validate),'sliced_pinned_validation','exec'),env)
    return env['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
