#include <metal_stdlib>
using namespace metal;

struct Point { float r,m,edge,b,eps,h,k,A,B; uint tuned; };
struct Row { Point point; float exact,leading,gap,offset; uint index,salt; };

inline uint hash32(uint x) {
    x ^= x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x,uint k) {
    return float(hash32(x+0x9e3779b9u*k)>>8)*0x1p-24f;
}
inline float logs(uint x,uint k,float a,float b) { return exp2(mix(a,b,unit(x,k))); }

inline float logone(float x) {
    if (abs(x)>0.03f) return log(1.0f+x);
    return x*(1.0f+x*(-0.5f+x*(1.0f/3.0f+x*(-0.25f+x*(0.2f-x/6.0f)))));
}
inline float expminus(float x) {
    if (abs(x)>0.08f) return exp(x)-1.0f;
    return x*(1.0f+x*(0.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}
inline float atanh_stable(float x) {
    const float xx=x*x;
    return xx<0.0025f ? x*(1.0f+xx*(1.0f/3.0f+xx*(0.2f+xx/7.0f)))
                     : 0.5f*log((1.0f+x)/(1.0f-x));
}

inline float small_sinh(float x) {
    if (abs(x) > 0.25f) return sinh(x);
    const float x2=x*x;
    return x*(1.0f+x2*(1.0f/6.0f+x2*(1.0f/120.0f+x2*(1.0f/5040.0f+x2/362880.0f))));
}

inline float3 panel(Point p,constant float2* nodes,uint n,float lo,float hi,bool evaluate) {
    if (hi<=lo) return float3(0);
    const float width=hi-lo;
    float3 total=0;
    for (uint j=0;j<n;++j) {
        const float t0=nodes[j].x;
        const float u=lo+width*t0*t0*t0;
        const float t=u/p.r,t2=t*t;
        const float c=p.edge*t2*t2;
        const float v=p.eps+c*(1.0f+p.b*t2);
        const float z=u*u/v;
        if (z>192.0f) continue;
        const float common=width*nodes[j].y*3.0f*t0*t0*rsqrt(v)*exp(-0.5f*(u*u+z));
        const float slope=c*t/v;
        if (evaluate) {
            const float x=p.h*slope, xx=x*x;
            const float le=-0.25f*logone(-xx)-0.5f*z*xx/(1.0f-xx);
            const float odd=0.5f*(z*x/(1.0f-xx)-atanh_stable(x));
            const float au=p.k*p.h*t;
            const float sh=small_sinh(0.5f*odd)/p.h;
            const float E=exp(le);
            const float residual=cosh(au)*(expminus(le)/(p.h*p.h)+2.0f*E*sh*sh)
                                -(small_sinh(au)/p.h)*E*(small_sinh(odd)/p.h);
            total.z += common*residual;
        } else {
            total.x += common*(z*z-6.0f*z+3.0f)*slope*slope/8.0f;
            total.y += common*t*(z-1.0f)*slope/2.0f;
        }
    }
    return total;
}

inline float3 integrate(Point p,constant float2* nodes,uint n,bool eval) {
    const float a=min(p.r,8.0f*sqrt(p.eps));
    const float b=max(a,min(p.r,0.25f*p.r*p.r/sqrt(p.edge)));
    const float c=max(b,min(p.r,p.r*p.r/sqrt(p.edge)));
    return panel(p,nodes,n,0,a,eval)+panel(p,nodes,n,a,b,eval)
         + panel(p,nodes,n,b,c,eval)+panel(p,nodes,n,c,p.r,eval);
}

kernel void quartic(constant float2* nodes [[buffer(0)]],device Row* output [[buffer(1)]],
                    constant uint4& params [[buffer(2)]],
                    device const uint* indices [[buffer(3)]],
                    uint id [[thread_position_in_grid]]) {
    const uint i=indices[id];
    const uint seed=hash32(i^0xd5bf242bu^params.w);
    const float r=logs(seed,1,-6.0f,0.0f);
    const float m=1.0f+logs(seed,2,-10.0f,3.0f);
    const float b=mix(-0.375f,0.75f,unit(seed,3));
    const float eps=logs(seed,4,-53.15084952f,-6.64385619f);
    const float edge=min(r*r*logs(seed,5,-5.0f,7.0f),0.8f*m/(1.18f+abs(b)));
    Point p=Point{r,m,edge,b,eps,0,0,0,0,0};
    const float3 ab=integrate(p,nodes,params.y,false);
    p.A=ab.x; p.B=ab.y;
    const float root=ab.x/ab.y;
    p.tuned=uint((i&1u)!=0 && abs(ab.y)>1e-4f && abs(root)<12.0f);
    const float dk=logs(seed,6,-12.0f,-2.0f);
    p.k=p.tuned!=0 ? root+((i&2u)!=0?dk:-dk) : mix(-12.0f,12.0f,unit(seed,7));
    p.h=min(logs(seed,8,-14.0f,-2.5f),3.0f*r/max(1.0f,abs(p.k)));
    const float minimum=integrate(p,nodes,params.y,true).z;
    const float a=p.k*p.h/r;
    const float target=m+eps;
    const float norm=2.0f*sqrt(target*0.1591549430918953f)*exp(-0.5f*a*a*(1.0f-1.0f/target));
    const float exact=norm*minimum, leading=norm*(ab.x-p.k*ab.y);
    const float metric=abs(exact-leading);
    output[id]=Row{p,exact,leading,metric,dk,i,params.w};
}
