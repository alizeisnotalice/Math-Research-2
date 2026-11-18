"""Certify each true first-hit class with a root bracket and one placement per root-containing angle."""
from fractions import Fraction as F
from flint import arb
from certify_ah_contact_slices import validate as test


def sliced_bin(ns,j,begin=0,end=75):
    candidates=[('14','23',300,1),('36','46',180,1),('34','46',180,1),('14','24',300,1),('14','25',300,1),('35','46',180,1),('36.5','46',180,1),('37','46',180,1),('14','24.5',300,1),('14','23.5',300,1),('35.5','46',180,1),('36.25','46',180,1),('36.75','46',180,1),('34.5','46',180,1)]
    rows=[]
    for ori in (0,1):
        for k in range(6):
            state=dict(ns);state['_first_orientation']=ori;state['_contact_levels']=(k,5-k)
            def block(lo,hi,preferred):
                opts=[preferred]+[v for v in candidates if v!=preferred]
                for l,h,t,s in opts:
                    mid=(lo+hi)//2
                    row=test(state,j,F(l),F(h),t,s,1,begin=mid,end=mid+1,angular=10,depth=2)
                    if not (arb(row['root_margin_lower'])>0 and arb(row['placement_margin_lower'])>0):continue
                    row=test(state,j,F(l),F(h),t,s,1,begin=lo,end=hi,angular=12,depth=3)
                    if arb(row['root_margin_lower'])>0 and arb(row['placement_margin_lower'])>0:
                        row.update(first_orientation=ori,contact_levels=[k,5-k]);rows.append(row)
                        return (l,h,t,s)
                if hi-lo<=1:
                    raise AssertionError((j,ori,k,lo,hi))
                mid=(lo+hi)//2;best=block(lo,mid,preferred);return block(mid,hi,best)
            start=candidates[0] if (ori and k!=5) or (not ori and k in (0,5)) else (candidates[2] if ori else candidates[1])
            block(begin,end,start)
            print('sliced',j,ori,k,len(rows),flush=True)
    return rows
