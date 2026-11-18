"""Expand every signed-area identity used by the rigid-translation screen in Z[x]."""
import json
from flint import fmpz_mpoly_ctx

def main():
    if not __debug__:
        raise SystemExit("Assertions required.")
    ring=fmpz_mpoly_ctx.get(tuple("ax ay bx by c d s t ux uy vx vy x y".split()))
    ax,ay,bx,by,c,d,s,t,ux,uy,vx,vy,x,y=ring.gens()
    a,b,u,v=(ax,ay),(bx,by),(ux,uy),(vx,vy)
    def det(p,q): return p[0]*q[1]-p[1]*q[0]
    def dot(p,q): return p[0]*q[0]+p[1]*q[1]
    def add(p,q): return (p[0]+q[0],p[1]+q[1])
    def move(p,co,si): return (co*p[0]-si*p[1],co*p[1]+si*p[0])
    ra,rb=move(a,c,s),move(b,d,t)
    lhs=det(add(u,ra),add(v,rb))
    expected=(det(u,v)+d*det(u,b)+t*dot(u,b)-c*det(v,a)-s*dot(v,a)
              +(c*d+s*t)*det(a,b)+(c*t-s*d)*dot(a,b))
    assert lhs==expected
    assert det(ra,move(b,c,s))==(c*c+s*s)*det(a,b)
    assert c*c+s*s==((c+1)**2+s*s-1)-2*c
    assert det(add(ra,u),(x,y))==ux*y-uy*x+c*(ax*y-ay*x)-s*(ax*x+ay*y)
    print(json.dumps(dict(distinct_body_edge=True,same_body_edge=True,same_body_circle_reduction=True,
                          first_variation=True,cosine_convention="cos(delta)-1",sine_convention="sin(delta)")))

if __name__=="__main__":
    main()
