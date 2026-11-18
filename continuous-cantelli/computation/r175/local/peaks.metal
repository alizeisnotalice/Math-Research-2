#include <metal_stdlib>
using namespace metal;

struct Shape { float c,L,s,w1,w2,h1,h2,l1,r1,l2,r2; };
struct Row {
    Shape shape;
    float D0,D1,D2,parab,det,sine,score,maxValue;
    uint mode,index,tested,witness,bad,salt;
};
constant float norm=.3989422804014327f;
constant float pi=3.1415926535897932f;
inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float pick(uint seed,uint i) {return float(hash32(seed+i*0x9e3779b9u)>>8)*0x1p-24f;}
inline float density(float x,float v) {return norm*rsqrt(v)*exp(-.5f*x*x/v);}
inline float4 rhs(Shape p,constant float2* nodes,uint order) {
    float4 out=0.f;
    for(uint j=0;j<order;++j) {
        const float u=nodes[j].x,x=p.c+(u-.5f)*p.L;
        const float w=nodes[j].y*density(x,p.s);
        out+=w*float4(1.f,u,u*u,sin(pi*u));
    }
    return out;
}
inline float4 peak(float c,float L,float left,float width,float lh,float rh,float h,
                   constant float2* nodes,uint order) {
    float4 out=0.f;
    const float start[3]={0.f,lh,rh},length[3]={lh,rh-lh,1.f-rh};
    for(uint i=0;i<3u;++i) {
        for(uint j=0;j<order;++j) {
            const float t=nodes[j].x, u=left+width*(start[i]+length[i]*t);
            const float x=c+(u-.5f)*L,b=h*(i==0u?t:(i==1u?1.0f:1.0f-t));
            const float w=width*length[i]*nodes[j].y*density(x,1.f);
            out+=w*float4(1.f,u,u*u-b,sin(pi*u)*exp(.5f*pi*pi*b));
        }
    }
    return out;
}
kernel void peaks_scan(constant float2* nodes [[buffer(0)]],device Row* out [[buffer(1)]],
                       constant uint4& pars [[buffer(2)]],uint id [[thread_position_in_grid]],
                       uint lane [[thread_index_in_simdgroup]]) {
    const uint index=pars.z+id,g=id/32u,mode=g%6u;
    const uint seed=hash32(index^0x178b2badu^pars.w);
    Shape p={mix(-1.2f,1.2f,pick(seed,1)),mix(.35f,2.4f,pick(seed,2)),
             mix(1.025f,2.05f,pick(seed,3)),mix(.03f,.38f,pick(seed,4)),
             mix(.03f,.38f,pick(seed,5)),mix(.06f,1.45f,pick(seed,6)),
             mix(.06f,1.45f,pick(seed,7)),mix(.045f,.38f,pick(seed,8)),
             mix(.62f,.955f,pick(seed,9)),mix(.045f,.38f,pick(seed,10)),
             mix(.62f,.955f,pick(seed,11))};
    if(mode==1u) {
        p.s=mix(1.025f,1.48f,pick(seed,3));
        p.h1=min(p.h1,(2.0f-p.s)/(p.L*p.L));
        p.h2=min(p.h2,(2.0f-p.s)/(p.L*p.L));
    }
    if(mode==2u) {
        p.c=0.0f;
        p.w2=p.w1;p.h2=p.h1;p.l2=1.0f-p.r1;p.r2=1.0f-p.l1;
    }
    if(mode==3u) {
        p.c=mix(-.16f,.16f,pick(seed,1));p.L=mix(.85f,1.15f,pick(seed,2));
        p.s=mix(1.14f,1.37f,pick(seed,3));
        p.w1=mix(.065f,.13f,pick(seed,4));p.w2=mix(.12f,.19f,pick(seed,5));
        p.h1=mix(.58f,.91f,pick(seed,6));p.h2=mix(.35f,.62f,pick(seed,7));
    }
    if(mode==4u) {
        p.w1=mix(.014f,.095f,pick(seed,4));p.w2=mix(.014f,.095f,pick(seed,5));
        p.h1=mix(.82f,2.4f,pick(seed,6));p.h2=mix(.82f,2.4f,pick(seed,7));
    }
    if(mode==5u) {
        p.c=0.0f;p.L=1.0f;p.s=1.125f+.03125f*float((seed>>27)%9u);
        p.w1=.0625f+.03125f*float((seed>>1)%5u);
        p.w2=.0625f+.03125f*float((seed>>4)%5u);
        p.h1=.375f+.0625f*float((seed>>9)%10u);
        p.h2=.375f+.0625f*float((seed>>14)%10u);
        p.l1=.125f;p.r1=.875f;p.l2=.125f;p.r2=.75f;
        if(lane==0u) p=Shape{0.f,1.f,1.25f,.09375f,.15625f,.75f,.5f,.125f,.875f,.125f,.75f};
    }
    const float4 target=rhs(p,nodes,pars.y);
    const float4 local=peak(p.c,p.L,0.f,p.w1,p.l1,p.r1,p.h1,nodes,pars.y)
                      +peak(p.c,p.L,1.f-p.w2,p.w2,p.l2,p.r2,p.h2,nodes,pars.y);
    const float4 D=target-local;
    const float parab=D.y-D.z,det=D.x*D.z-D.y*D.y;
    const float cap=p.s+p.L*p.L*max(p.h1,p.h2);
    const bool good=all(isfinite(D))&&isfinite(det);
    const bool witness=good&&D.x>0.f&&det>0.f&&parab>0.f&&D.w<0.f;
    const float score=good?min(min(D.x,parab),min(det/max(D.x,0.001f),-D.w)):-INFINITY;
    const float top=simd_max(score);
    const uint winner=simd_min(score==top?lane:32u);
    const uint count=simd_sum(uint(witness)),bad=simd_sum(uint(!good));
    if(lane==winner) out[g]=Row{p,D.x,D.y,D.z,parab,det,D.w,score,cap,mode,index,32u,count,bad,pars.w};
}
