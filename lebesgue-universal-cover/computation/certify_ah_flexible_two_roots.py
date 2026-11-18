"""Select a support-sign class and allow a free multiple of the root centering equation in any edge margin."""
from fractions import Fraction as F
from flint import arb
from certify_ah_affine_root_class import validate as conditional


def either_root(ns,j):
    options=[(F('23.4'),F('35.4'),F('.00001')),(F('23.5'),F('35.5'),F('.00001')),(F('23.3'),F('35.3'),F('.000005')),(F('23.4'),F('35.6'),F('.000005')),(F(24),F(36),F('.00004'))]
    rows=[]
    def check(lo,hi,params):
        anchor,other,tau=params
        tested=[]
        for l,r,t,sgn in ((F(14),anchor,300,-1),(other,F(46),180,1)):
            row=conditional(ns,j,l,r,t,1,1,begin=lo,end=hi,angular=12,depth=4,anchor=anchor,threshold=tau,condition_sign=sgn)
            if arb(row['root_margin_lower'])<=0 or arb(row['placement_margin_lower'])<=0:return None
            tested.append(row)
        return tested
    def block(lo,hi,first):
        for opt in [first]+[item for item in options if item!=first]:
            if hi-lo>2 and not check((lo+hi)//2,(lo+hi)//2+1,opt):continue
            got=check(lo,hi,opt)
            if got:
                rows.extend(got);return opt
        if hi-lo==1:raise AssertionError((j,lo,hi,'split support classes'))
        mid=(lo+hi)//2;first=block(lo,mid,first);return block(mid,hi,first)
    block(0,75,options[0])
    return {'reference_bin':j,'two_bracket_classes':rows,'root_margin_lower':str(min(arb(q['root_margin_lower']).lower() for q in rows)),'placement_margin_lower':str(min(arb(q['placement_margin_lower']).lower() for q in rows)),'placed_leaves':sum(q['placed_leaves'] for q in rows),'root_free_leaves':sum(q['root_free_leaves'] for q in rows)}
