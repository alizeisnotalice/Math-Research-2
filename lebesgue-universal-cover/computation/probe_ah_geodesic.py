"""Numerical test of the L(theta) boundary arc, not a universal-cover proof."""
import json,sys
from math import acos,atan2,sqrt,pi
from reconstruct_ah import definition8,add,sub,mul,norm,dot,cross,rot
from reconstruct_s_primitives import circle_line_intersections,arc_green_integral

def run(sigma,k):
    a=definition8(sigma);th=a['theta'];tau=a['tau'];Y=a['Y'];Z=a['Z'];C=a['L'](th)
    z=sub(Z,C);d2=dot(z,z);assert d2>1
    tangent=[add(C,add(mul(1/d2,z),mul(sign*sqrt(d2-1)/d2,(-z[1],z[0])))) for sign in (-1,1)]
    T=min(tangent,key=lambda t:norm(sub(t,Y)))
    u=sub(Y,C);v=sub(T,C);angle=atan2(cross(u,v),dot(u,v))
    path=[add(C,rot(u,angle*j/k)) for j in range(k+1)]
    path += [add(T,mul(j/k,sub(Z,T))) for j in range(1,k+1)]
    def inverse(x):
        L=min(circle_line_intersections(x,(.5,-sqrt(3)/2)),key=lambda l:norm(sub(l,a['L'](0))))
        N=min(circle_line_intersections(x,(-.5,-sqrt(3)/2)),key=lambda n:norm(sub(n,a['M'])))
        s=acos(max(-1,min(1,L[0]+.5)))-pi/6
        t=dot(sub(N,a['M']),a['e_n'])
        return s,t
    st=[inverse(x) for x in path]
    residual=max(norm(sub(x,a['X'](s,t))) for x,(s,t) in zip(path,st))
    overshoot=max(max(abs(s)-th,abs(t)-tau,0) for s,t in st)
    paired=[a['X'](-s,-t) for s,t in st]
    chain=path+list(reversed(paired))[1:]
    green=sum(cross(x,y) for x,y in zip(chain,chain[1:]))
    green += arc_green_integral(a['P'],a['X_as'],a['q']['G'])
    green += arc_green_integral(a['X_as'],Y,a['d']['E3'])
    result={'sigma_degrees':sigma,'segments_per_piece':k,'Y_distance_G':norm(sub(Y,a['q']['G'])),'Y_distance_Ltheta':norm(sub(Y,C)),'Z_distance_Ltheta':sqrt(d2),'arc_angle':angle,'tangent':T,'tangent_dot_residual':abs(dot(sub(T,C),sub(Z,T))),'inverse_residual':residual,'parameter_box_overshoot':overshoot,'paired_endpoint_residual':max(norm(sub(paired[0],a['P'])),norm(sub(paired[-1],Z))),'AH_area_numerical':abs(green/2)}
    assert residual<1e-11 and overshoot<1e-11
    assert result['tangent_dot_residual']<1e-11 and result['paired_endpoint_residual']<1e-11
    return result
if __name__=='__main__':
    print(json.dumps({'python':sys.version,'scope':'sampled coordinate/geodesic candidate only; no interval or universal-cover certificate','runs':[run(1.5494,k) for k in (100,1000,4000)]},indent=2))
