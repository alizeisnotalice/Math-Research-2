"""Normalize by the actual larger graph deficit and track the first hit's straight/paired label."""
import contextlib,inspect,io
from fractions import Fraction as F
from certify_ah_margin_envelope import replace_line
from certify_ah_centered_deficit import geometry as unlabelled_geometry
from certify_ah_single_contact_chain import support_bound as unlabelled_support,local_slope
from certify_ah_attainable_root import validate as root_validation


def geometry(archive,extra,*,far,negative,eps):
    template=inspect.getsource(unlabelled_geometry)
    marker="    ns={'__file__':str(archive/'certify_ah_two_line_cap.py'),'__name__':'__main__'}"
    prefix="""    prefix=replace_line(prefix,'assert hmax<A("0.000080")','assert hmax<A("0.00016")')
    prefix=replace_line(prefix,'parameter_boxes=[]','labelled_deficits=[[zero for _ in ref_boxes] for _ in range(2)]\\nparameter_boxes=[]')
    prefix=replace_line(prefix,'    parameter_boxes.append((l,pl,used,a))','''    for orientation,ee in enumerate((l,pl)):
        if not ee>0:
            for k,rr in enumerate(ref_boxes):
                if ee.overlaps(rr):labelled_deficits[orientation][k]=max(labelled_deficits[orientation][k],cell_deficit)
    parameter_boxes.append((l,pl,used,a))''')
"""
    template=replace_line(template,marker,prefix+marker);sc=dict(unlabelled_geometry.__globals__);exec(compile(template,'first_hit_labelled_geometry','exec'),sc);ns,code,old=sc['geometry'](archive,extra,far=far,negative=negative,eps=eps);assert all(ns['deficit_by_bin'][j]==max(row[j] for row in ns['labelled_deficits']) for j in range(len(ns['ref_boxes'])));return ns,code,old


def support_bound(ns,j,a,b):
    template=inspect.getsource(unlabelled_support)
    mark='    # Only the coupled bound is exported.'
    extra='''    src=replace_line(src,"first+de-partner+support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de,-first),zero,A(2))", "-first+de-partner+support_chain((x/(m*A('.95')),-eta/m,-xi/m,partner,-de),zero,A(2))")
    src=replace_line(src,'return extremal((value(ll,lr,a2,t4,dt,t2,local_m[0]),value(lr,ll,a4,t2,-dt,t4,local_m[1])))',"return extremal((zero,ns['labelled_deficits'][0][j]/E*value(ll,lr,a2,t4,dt,t2,local_m[0]),ns['labelled_deficits'][1][j]/E*value(lr,ll,a4,t2,-dt,t4,local_m[1])))")
'''
    template=replace_line(template,mark,extra+mark);sc=dict(unlabelled_support.__globals__);exec(compile(template,'first_hit_simplex','exec'),sc);return sc['support_bound'](ns,j,a,b)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(root_validation);sc={**root_validation.__globals__,'support_bound':support_bound};exec(compile(src,'actual_first_hit_root','exec'),sc);return sc['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
