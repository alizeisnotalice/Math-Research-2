"""Direct finite Green formula for the forced and rejected contours."""
from fractions import Fraction as F

def geometry(s, arc_checks=None):
    from flint import arb_series
    from certify_low_sigma_ah_budget import model,b,g,sub,add,dot,rot,cross,mul
    from certify_three_sigma_area_budget import exact
    n=b.directions
    m=[g.unit(g.PI/6+s+i*g.PI/3) for i in range(6)]
    meet=g.line_pair
    half=g.HALF
    def proj(c,v):return sub(c,mul(dot(c,v)-half,v))
    b1,f1,d1=[meet(n[i],n[j]) for i,j in ((0,1),(2,3),(4,5))]
    B,C,X,Y0=[meet(n[i],m[j]) for i,j in ((1,0),(0,5),(4,3),(3,2))]
    G=g.circle_line(Y0,n[5],1)
    r=meet(m[2],n[1],half,0)
    S=g.circle_line(r,n[0],-1)
    T=g.circle_pair(S,B,1)
    E=sub(S,n[0])
    Q,P,I,K=[proj(a,v) for a,v in ((X,n[1]),(G,n[2]),(B,m[3]),(Y0,m[5]))]
    c=(-X[0],X[1]);yy=(c[0]*0,c[1]+(1-c[0]**2).sqrt())
    D=model.D_of(s);ta=add(b.z,D)
    cp=add(b.z,mul(exact(F('.435')),D))
    ws=[F('.435')]+[F(i,20) for i in range(9,21)]
    qs=[b.anchors(D,exact(w))[6] for w in ws]
    def arc(a,z,c):
        aa,vv=sub(a,c),sub(z,c)
        cr,dt=cross(aa,vv),dot(aa,vv)
        zeroth=lambda v:v[0] if isinstance(v,arb_series) else v
        assert zeroth(dt)>0
        if arc_checks is not None: arc_checks.append(cr)
        return cross(c,sub(z,a))+(cr/dt).atan()
    chain=[ta,cp]+qs
    base=(cross(b1,Q)+arc(Q,yy,X)+arc(yy,ta,c)
          +sum(cross(a,z) for a,z in zip(chain,chain[1:]))
          +arc(qs[-1],P,G)+cross(P,f1)+cross(f1,E)
          +arc(T,I,B)+cross(I,X)+cross(X,d1)+cross(d1,G)+arc(G,K,Y0)+cross(K,C)+cross(C,b1))/2
    return base+arc(E,T,S)/2, base+cross(E,T)/2
