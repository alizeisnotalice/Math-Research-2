"""Certify either of two disjoint angle brackets selected by a sign of the same body's support."""
from fractions import Fraction as F
from flint import arb
from certify_ah_conditional_second_root import validate as conditional
from certify_ah_local_physical_support import validate as ordinary,support_bound


def either_root(ns,j):
    A=ns['A'];anchor=F(23);threshold=F('.00003');start=F(14)
    first=ordinary(ns,j,start,anchor,300,1,1,angular=12,depth=4)
    assert arb(first['placement_margin_lower'])>0,first
    n0=ns['av'](ns['d']['n'][0]);terms=[(A((-1)**k),ns['rotate'](n0,ns['trig'](A(start+60*k)))) for k in range(3)]
    left=A(1);joint=A(1);ri=None;ji=None
    six=[(A(s*(-1)**k),ns['rotate'](n0,ns['trig'](A(h+60*k)))) for s,h in ((1,anchor),(-1,F(36))) for k in range(3)]
    for i in range(75):
        r0,r1=ns['knots'][j:j+2];step=(r1-r0)/75;a=r0+i*step;b=a+step
        E,v,e,bound,_=support_bound(ns,j,a,b)
        margin=(ns['F_at'](v,e,start)-E*bound([(-c,u) for c,u in terms])).lower()
        diff=(ns['F_at'](v,e,F(36))-ns['F_at'](v,e,anchor)-E*bound(six)).lower()
        if margin<left:left,ri=margin,i
        if diff<joint:joint,ji=diff,i
    assert left>0 and joint>A(threshold),(j,left,joint)
    first['root_margin_lower']=str(min(left,A(threshold)));first['root_id']=[ri,str(start)]
    first['anchor_upper_condition']=str(-threshold);first['left_endpoint_unconditional_margin']=str(left)
    second=conditional(ns,j,F(36),F(46),180,1,1,angular=12,depth=4,anchor=anchor,threshold=threshold)
    assert arb(second['root_margin_lower'])>0 and arb(second['placement_margin_lower'])>0,second
    rows=(first,second)
    return {'reference_bin':j,'two_bracket_classes':rows,'F36_minus_F23_lower':str(joint),'F36_minus_F23_min_index':ji,'root_margin_lower':str(min(arb(q['root_margin_lower']).lower() for q in rows)),'placement_margin_lower':str(min(arb(q['placement_margin_lower']).lower() for q in rows)),'placed_leaves':sum(q['placed_leaves'] for q in rows),'root_free_leaves':sum(q['root_free_leaves'] for q in rows)}
