"""Use the three genuine contact pairs of the initial hexagon to enclose the true compatibility equation."""
import inspect
from certify_ah_attainable_root import validate as base_validate
from certify_ah_root_pose_choices import validate as pose_validate
from certify_ah_margin_envelope import replace_line


def pinned_slope(ns,v,E):
    A=ns['A']
    normals=list(map(ns['av'],ns['d']['n']))
    # v1,v3,v0 belong to the sides with outward normals -n0,-n1,-n2, respectively.
    tangent_sum=sum((-1)**k*ns['dot'](v[i],(-normals[k][1],normals[k][0])) for k,i in enumerate((1,3,0)))
    # Two L contacts move <= (2/3)/.95 E and (2/3) E, in either order, and the N contact <=E.
    radius=(A(2)/(3*A('.95'))+A(2)/3+1)*E
    return tangent_sum+(-radius).union(radius)


def pinned_equation(ns,slope,angle):
    c,s=ns['trig'](angle)
    # Use 1-cos(t)=2*sin(t/2)^2 to avoid losing a sign in small intervals.
    _,half_sin=ns['trig'](angle/2)
    radius=3*half_sin**2
    return slope*s+(-radius).union(radius)


def validate(ns,j,left,right,turn,sign,orientation,begin=0,end=75,angular=12,depth=5):
    src=inspect.getsource(base_validate)
    src=replace_line(src,"        def root_terms(angle):","        slope=pinned_slope(ns,v,E)\n        def root_terms(angle):")
    src=replace_line(src,"            if slack<rootmin:","            slack=max(slack,(s*pinned_equation(ns,slope,A(phi))).lower())\n            if slack<rootmin:")
    src=replace_line(src,"            if no_root>0:","            pin=pinned_equation(ns,slope,interval);no_root=max(no_root,pin.lower(),(-pin).lower())\n            if no_root>0:")
    template=inspect.getsource(pose_validate)
    template=replace_line(template,"src=inspect.getsource(base_validate)","src=pinned_source")
    template=replace_line(template,"'support_bound':support_bound","'support_bound':support_bound,'pinned_slope':pinned_slope,'pinned_equation':pinned_equation")
    sc={**pose_validate.__globals__,'pinned_source':src,'pinned_slope':pinned_slope,'pinned_equation':pinned_equation}
    exec(compile(template,'pinned_dihedral_root','exec'),sc)
    return sc['validate'](ns,j,left,right,turn,sign,orientation,begin,end,angular,depth)
