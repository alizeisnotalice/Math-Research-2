"""Uniform independent checks for the AH fourth support and the mixed-area term."""
from fractions import Fraction as F
import inspect
from flint import arb, arb_series

def strengthened_root():
    from certify_sigma_paired_rigidity import rig
    source=inspect.getsource(rig)
    old='if gap>0:'
    assert source.count(old)==1
    # The returned new interior root margin is for D alone, after demanding
    # enough room against the exact amplitude bound for the Reuleaux triangle.
    source=source.replace(old,'if gap>g.exact(F(24,100))*g.exact(F(4,100000))/(1-g.exact(F(4,100000))):')
    ns=rig.__globals__.copy()
    exec(compile(source,'rig_with_triangle_drift','exec'),ns)
    return ns['rig']

def fourth(l,r,sig,rig):
    from certify_low_sigma_ah_budget import model,exact,b,g,add,sub,mul,dot,rot,cross
    w=F('.435');tau=exact(F(4,100000))
    n,stats=rig(l,r,w,w,endcap=F('.1'),min_cell=F(1,10**7))
    assert (1-tau)*stats[3]-3*tau>0 and (1-tau)*stats[4]-exact('.24')*tau>0
    def replacement(s):
        D=model.D_of(s);ta=add(b.z,D)
        x=add(b.z,mul(exact(w),D));q=b.anchors(D,exact(w))[6]
        upper_left=(-g.HALF,1/(2*arb(3).sqrt()))
        xt=add(mul(1-tau,x),mul(tau,upper_left))
        return (cross(ta,xt),cross(xt,q),cross(sub(ta,q),sub(x,upper_left))/2)
    mid=b.degree((l+r)/2)
    aa,bb,cost=[v+ds[1]*(sig-mid) for v,ds in zip(replacement(mid),replacement(arb_series([sig,1],2)))]
    assert min(aa,bb)>0 and 0<cost<exact('.00192'),(l,r,aa,bb,cost)
    def tip_slack(a,z):
        de=model.D_of(b.degree(a).union(b.degree(z)))
        qu=b.anchors(de,arb_series([exact(w),1],2))[6]
        vel=tuple(v[1] for v in qu);p=tuple(v[0] for v in qu)
        val=-dot(rot(vel),sub((arb(0),1/arb(3).sqrt()),p))/g.norm(vel)
        if val>exact('.009'):return val.lower()
        assert z-a>F(1,10**8),(a,z,val)
        mm=(a+z)/2
        return min(tip_slack(a,mm),tip_slack(mm,z))
    return n,cost.upper(),tip_slack(l,r)

def mixed_area(l,r):
    from certify_low_sigma_ah_budget import model,exact,b,g,add,sub,mul,dot
    from certify_sigma_mandatory_arcs import anchors
    g.SIGMA=b.degree(l).union(b.degree(r));n,m,*_=g.setup()
    ti=anchors('TI')[0];gk=anchors('GK')[0];qy=anchors('QY')[0]
    tips=[g.line_pair(n[i],n[j]) for i,j in ((0,1),(2,3),(4,5))]
    c=g.line_pair(n[0],m[5])
    points=tips+[c,(-c[0],c[1]),qy['C'],gk['G'],ti['T'],ti['E'],qy['Q'],
        sub(gk['G'],n[5]),gk['K'],qy['Y'],add(b.z,model.D_of(g.SIGMA))]
    # G-n5=P lies on the T'P disk arc (G-n5 = projection of G on n2).
    # I = project(B3,m3) (the TI endpoint).
    p=ti['B'];points.append(sub(p,mul(dot(p,m[3])-g.HALF,m[3])))
    total=arb(0);count=4096
    for j in range(count):
        a=2*g.PI*j/count;z=2*g.PI*(j+1)/count
        theta=a.union(z)
        u=g.unit(theta);h=max(g.HALF,max(dot(u,z).lower() for z in points))
        mass=(z-a)/2+((3*z).cos()-(3*a).cos())/6
        total+=h*max(arb(0),mass.lower())
    v=total/2
    assert v>exact('.8223'),(l,r,v)
    return v.lower()
