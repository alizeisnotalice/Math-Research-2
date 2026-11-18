"""Disk and one smooth width-one body exclude octagons with eight moving normals."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import json,sys
from flint import arb,ctx,__version__ as flint_version
ctx.prec=160
def q(x):
    x=F(x);return arb(x.numerator)/x.denominator
def deg(x):return q(x)*arb.pi()/180
def unit(t):return (t.cos(),t.sin())
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def norm(a):return (a[0]**2+a[1]**2).sqrt()
def vector_sum(X,w):return [sum(v[i]*y for v,y in zip(X,w)) for i in range(2)]
def bits(x):
    m,e=x.man_exp()
    return str(F(int(m))*(F(2)**int(e)))
def run(delta_min,delta_max,normal_error,panels):
    stats={}; rows=[]
    B=q('.844082046521218581802948589770767704354')
    h=arb(1)/2;blend=q('.002')
    patterns=[((1,3),(0,1,1,0,0,1)),((1,3),(0,0,1,0,0,1)),((1,4),(0,0,1,0,0,1)),
              ((0,6),(0,0,1,0,0,0)),((5,6),(1,0,0,1,0,0)),((6,7),(1,0,0,1,0,1))]
    def rec(name,value):
        lo=value.lower();stats[name]=min(stats.get(name,lo),lo)
    for panel in range(panels):
        a=delta_min+(delta_max-delta_min)*panel/panels;b=delta_min+(delta_max-delta_min)*(panel+1)/panels
        d=deg(a).union(deg(b))
        er=(-deg(normal_error)).union(deg(normal_error))
        t=[o+er for o in (deg(0),deg(64),deg(96),deg(128)-d,deg(180),deg(212),deg(244),deg(308)+d)]
        N=list(map(unit,t));Z=[unit(5*x) for x in t];gaps=[v-u for u,v in zip(t,t[1:]+[t[0]+2*arb.pi()])]
        L=[h*((gaps[i-1]/2).tan()+(gaps[i]/2).tan()) for i in range(8)]
        A=sum(L)/4
        assert min(L)>q('.20') and min(gaps)>0 and max(gaps)<arb.pi()/2
        # Dimensionless bounded equilibria: start with six fixed dual coordinates,
        # solve the normal balance for the remaining two, then move slightly inside
        # the valid box via the facet-length equilibrium.
        first=[]
        for (i,j),mask in patterns:
            f=[k for k in range(8) if k not in (i,j)]
            lam=[arb(0)]*8
            for k,bit in zip(f,mask):lam[k]=L[k]*bit
            rhs=[-x for x in vector_sum(N,lam)]
            det=cross(N[i],N[j])
            lam[i]=cross(rhs,N[j])/det;lam[j]=cross(N[i],rhs)/det
            w=[(1-blend)*v+blend*ell/2 for v,ell in zip(lam,L)]
            assert min(w)>0 and min(ell-u for ell,u in zip(L,w))>0,('dual',panel,i,j)
            first.append(w)
        dual=first+[[ell-x for ell,x in zip(L,w)] for w in first]
        V=[[q('999/48000')*z for z in vector_sum(Z,w)] for w in dual]
        radii=[cross(v,z)/norm(sub(z,v)) for v,z in zip(V,V[1:]+V[:1])]
        r=min(x.lower() for x in radii)
        assert r>0,('radius',panel,r)
        simple=A+r/2
        # For X = area(P)-A, the triangle on facet i bounds its actual
        # support surplus by 2X/L_i. Positive adjacent edge lengths persist
        # for parallel interpolation whenever 2*gamma*X < 1.
        gamma=max(((1/gaps[i-1].tan()+1/gaps[i].tan())/L[i]**2).upper() for i in range(8))
        assert gamma*r<1
        wall=A+r/(1+gamma*r)
        assert wall>B,(panel,a,b,A,r,gamma,wall)
        rec('base_area',A);rec('negative_base_area_upper',-A)
        rec('inradius',r);rec('arbitrary_extension_lower',simple)
        rec('octagon_area',wall);rec('B_excess',wall-B);rec('inverse_edge_condition_lower',1/gamma)
        # Unperturbed baseline alone, not the independently angled intervals:
        # an upward vertex at 34 degrees, then centered in exact strips 0,64.
        RT=arb(3).sqrt();h0=1-deg(26).cos()/RT
        bt=[deg(96),deg(128)-d,deg(212),deg(308)+d]
        hs=[1-deg(2).cos()/RT,(deg(26)+d).cos()/RT,1-deg(2).cos()/RT,1-(deg(26)-d).cos()/RT]
        slack=min((h-ah-(deg(64)-phi).sin()/deg(64).sin()*(h-h0)).lower() for phi,ah in zip(bt,hs))
        assert slack>0;rec('unperturbed_RT_cap_slack',slack)
        rows.append([str(a),str(b),bits(r.lower()),bits((wall-B).lower())])
    return dict(complete=True,d_degrees=[str(delta_min),str(delta_max)],independent_normal_errors_deg=str(normal_error),
        panels=panels,base=[0,64,96,'128-d',180,212,244,'308+d'],harmonic=5,amplitude='999/48000',B=str(F('.844082046521218581802948589770767704354')),
        lower_bounds={k:bits(v) for k,v in stats.items()},balanced_vector_patterns=patterns,panel_bounds=rows,python=sys.version,flint=flint_version,
        sha256=sha256(Path(__file__).read_bytes()).hexdigest(),precision=ctx.prec)

if __name__=='__main__':
    a,b,eps=map(F,sys.argv[1:4]);v=run(a,b,eps,int(sys.argv[4]))
    Path(sys.argv[5]).write_text(json.dumps(v,indent=1)+'\n');print({k:float(F(x)) for k,x in v['lower_bounds'].items()})
