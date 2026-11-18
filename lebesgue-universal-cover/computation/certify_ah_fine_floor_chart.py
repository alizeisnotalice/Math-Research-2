"""Use the true lower arc together with the finely subdivided contact tube and the cell's own first-hit budget."""
from certify_ah_physical_floor import FLOOR
from certify_ah_fine_contact_chart import fine_chart
from certify_ah_two_hit_budgets import geometry as complete
from certify_ah_margin_envelope import replace_line
import inspect


def floor_chart(archive,extra,*,far,negative,eps):
    old,src,_=fine_chart(archive,extra,far=far,negative=negative,eps=eps)
    src=replace_line(src,'parameter_boxes=[]',FLOOR+'parameter_boxes=[]')
    before='''    own_height = height_bound(a)
    x=add(add(z,mul(w,D)),mul(box(0,own_height.upper()),normal))
    cl,ct=chart(add(z,mul(A((a+b)/2),D)))'''
    after='''    if A(-b)>physical_rights[-1]:
        continue
    own_height=height_bound(a)
    floor=physical_floor(a,b)
    if own_height<floor:
        continue
    x=add(add(z,mul(w,D)),mul(floor.union(own_height.upper()),normal))
    cl,ct=chart(add(add(z,mul(A((a+b)/2),D)),mul(floor,normal)))'''
    src=replace_line(src,before,after)
    src=replace_line(src,'ht=box(0,own_height.upper())','ht=box(0,(own_height-floor).upper())')
    return old,src,None


def geometry(archive,extra,*,far,negative,eps):
    src=inspect.getsource(complete);scope={**complete.__globals__,'first_geometry':floor_chart}
    exec(compile(src,'fine_floor_first_hit_chart','exec'),scope)
    return scope['geometry'](archive,extra,far=far,negative=negative,eps=eps)
