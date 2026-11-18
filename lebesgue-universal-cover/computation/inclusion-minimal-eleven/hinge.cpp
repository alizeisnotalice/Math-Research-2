// Verify infeasibility of translations by integer hinge inequalities for two polygons.
#ifdef NDEBUG
#error This proof verifier requires assertions.
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using L=std::int64_t;
using W=__int128_t;
struct P{L x,y;};
using Poly=std::vector<P>;
struct Cut{L a,b,c;}; // accepts ax+by<=c in integer grid coordinates, translation real
struct Vertex{int i,j; W x,y,d;};
L det(P u,P v){return u.x*v.y-u.y*v.x;}
L area(const Poly& h){L a=0;for(size_t j=0;j<h.size();++j)a+=det(h[j],h[(j+1)%h.size()]);return a;}
Poly read_poly(std::istream& f){
    int n=0;f>>n;assert(n>=3 && n<10000);Poly h(n);
    for(auto& v:h){f>>v.x>>v.y;assert(f && std::abs(v.x)<=10000000 && std::abs(v.y)<=10000000);}
    for(int i=0;i<n;++i){P a=h[i],b=h[(i+1)%n],c=h[(i+2)%n];assert(det({b.x-a.x,b.y-a.y},{c.x-b.x,c.y-b.y})>0);}
    for(int i=0;i<n;++i){
        P a=h[i],b=h[(i+1)%n];
        for(int j=0;j<n;++j)if(j!=i && j!=(i+1)%n){
            P c=h[j];assert(det({b.x-a.x,b.y-a.y},{c.x-a.x,c.y-a.y})>0);
        }
    }
    assert(area(h)>0);return h;
}
Vertex meet(const std::vector<Cut>& cuts,int i,int j){
    auto [a,b,c]=cuts[i];auto [u,v,w]=cuts[j];W d=W(a)*v-W(u)*b;
    assert(d);W x=W(c)*v-W(w)*b,y=W(a)*w-W(u)*c;
    if(d<0){d=-d;x=-x;y=-y;}
    return {i,j,x,y,d};
}
bool inside(const Cut& c,const Vertex& v){return W(c.a)*v.x+W(c.b)*v.y<=W(c.c)*v.d;}
// A returned contradiction is certified anew from two or three input cuts.
// The polygon clipping is only a selector for these exact Farkas weights.
bool contradiction(const std::vector<Cut>& cuts){
    auto z=cuts.back();P nz{z.a,z.b};
    for(size_t i=0;i+1<cuts.size();++i){
        auto v=cuts[i];P nv{v.a,v.b};L u=det(nz,nv);
        if(!u && W(nz.x)*nv.x+W(nz.y)*nv.y<0){
            L w1=std::abs(nz.x?nz.x:nz.y),w2=std::abs(nz.x?nv.x:nv.y);
            if(W(w1)*v.c+W(w2)*z.c<0)return true;
        }
        for(size_t j=i+1;j+1<cuts.size();++j){
            auto w=cuts[j];P nw{w.a,w.b};L k=det(nv,nw),l=det(nw,nz);
            if((u>0 && l>0 && k>0) || (u<0 && l<0 && k<0))
                if(W(std::abs(l))*v.c+W(std::abs(u))*w.c+W(std::abs(k))*z.c<0)return true;
        }
    }
    return false;
}
struct Hinge{L a,b,c;};
void append(const Poly& self,const Poly& other,bool first,std::vector<Hinge>& z){
    for(size_t i=0;i<self.size();++i){
        auto p=self[i],q=self[(i+1)%self.size()];L a=q.y-p.y,b=p.x-q.x,v=a*p.x+b*p.y;
        L w=INT64_MIN;for(auto r:other)w=std::max(w,a*r.x+b*r.y);
        z.push_back({a,b,first?w-v:v-w});
    }
}
bool exclude(const Poly& h,const Poly& q,L budget){
    std::vector<Hinge> rows;append(h,q,true,rows);append(q,h,false,rows);
    // When the target is <=1 and the base contains a .48 disk, a point at
    // distance 3 would already create a triangle of area >1. The point polygon
    // contains the origin, hence t is inside [-3S,3S]^2.
    constexpr L T=30000000;
    std::vector<Cut> cuts={{1,0,T},{0,1,T},{-1,0,T},{0,-1,T}};
    std::vector<Vertex> v;for(int i=0;i<4;++i)v.push_back(meet(cuts,i,(i+1)%4));
    for(int steps=0;steps<65 && !v.empty();++steps){
        long double x=0,y=0;for(auto t:v){x+=(long double)t.x/t.d;y+=(long double)t.y/t.d;}x/=v.size();y/=v.size();
        Cut c{0,0,budget};
        for(auto [a,b,k]:rows)if(a*x+b*y+k>0){c.a+=a;c.b+=b;c.c-=k;}
        // Any selected subset is a valid affine minorant, independent of floating choice.
        if(!c.a && !c.b)return c.c<0;
        bool all=true;for(auto& t:v)all&=inside(c,t);
        if(all)return false;
        int id=(int)cuts.size();cuts.push_back(c);std::vector<Vertex> next;
        for(size_t j=0;j<v.size();++j){
            auto u=v[j],w=v[(j+1)%v.size()];bool I=inside(c,u),J=inside(c,w);
            if(I)next.push_back(u);
            if(I!=J){
                int edge=(u.i==w.i || u.i==w.j)?u.i:u.j;
                assert(edge==w.i || edge==w.j);
                next.push_back(meet(cuts,edge,id));
            }
        }
        if(next.size()==1 || next.size()==2)return false;
        v.swap(next);
    }
    return v.empty() && contradiction(cuts);
}
int main(int argc,char**argv){
    assert(argc==4 || argc==5);bool first_only=argc==5;std::ifstream host(argv[1]),bank(argv[2]);L S,T,target=std::stoll(argv[3]);int m,n;
    host>>S>>m;bank>>T>>n;assert(S==10000000 && T==S && m>=1 && n>=1 && n<=65536 && 0<target && target<=1000000);
    std::vector<Poly> q(2*n);for(auto& p:q){
        p=read_poly(bank);
        for(size_t j=0;j<p.size();++j)assert(det(p[j],p[(j+1)%p.size()])>0);
    }
    for(int i=0;i<m;++i){
        auto p=read_poly(host);
        for(size_t j=0;j<p.size();++j){
            P a=p[j],b=p[(j+1)%p.size()];L t=det(a,b),x=a.x-b.x,y=a.y-b.y;
            assert(t>0 && W(t)*t>W(4800000)*4800000*(x*x+y*y));
        }
        L e=2*(S*S/1000000)*target-area(p);
        int done=0,first=-1;
        if(e<0)done=2*n;
        else for(int j=0;j<2*n;++j){bool gr=exclude(p,q[j],e);if(gr)++done;else if(first<0){first=j;if(first_only)break;}}
        std::cout<<i<<" "<<done<<" "<<2*n<<" "<<first<<std::endl;
    }
}
