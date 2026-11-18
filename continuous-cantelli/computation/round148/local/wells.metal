#include <metal_stdlib>
using namespace metal;

struct Row {
    float eps,turns,sine,rawSine,first,second,stableScaled,rawScaled;
    float secondReal,secondImag,mass2,quadFloor,score;
    uint n,mode,band,tested,index,bad,salt;
};
constant float invSqrt2Pi=0.3989422804014327f;
constant float piValue=3.1415926535897932f;
constant float lead=0.12320154025412866f;
constant float secondCos=0.05005062572823977f;
inline uint hash32(uint x) {
    x^=x>>16; x*=0x7feb352du; x^=x>>15; x*=0x846ca68bu; return x^(x>>16);
}
inline float unit(uint x) {return float(hash32(x)>>8)*0x1p-24f;}
inline float log2rem(float x) {
    if (abs(x)>.12f) return (log(1.0f+x)-x)/(x*x);
    return -.5f+x*(1.0f/3.0f+x*(-.25f+x*(.2f+x*(-1.0f/6.0f+x*(1.0f/7.0f-x/8.0f)))));
}
inline float2 mul(float2 a,float2 b) {return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
inline float2 exp2norm(float2 fac,float eps) {
    const float2 z=fac*eps;
    if (max(abs(z.x),abs(z.y))>.25f) {
        const float r=exp(z.x);
        return (float2(r*cos(z.y)-1.0f,r*sin(z.y))-z)/(eps*eps);
    }
    const float2 c=float2(0.5f,0)+mul(z,float2(1.0f/6.0f,0)+mul(z,float2(1.0f/24.0f,0)
        +mul(z,float2(1.0f/120.0f,0)+mul(z,float2(1.0f/720.0f,0)+z/5040.0f))));
    return mul(mul(fac,fac),c);
}
inline float exp2real(float a,float eps) {
    const float z=a*eps;
    if (abs(z)>.25f) return (exp(z)-1.0f-z)/(eps*eps);
    return a*a*(.5f+z*(1.0f/6.0f+z*(1.0f/24.0f+z*(1.0f/120.0f+z*(1.0f/720.0f+z/5040.0f)))));
}

struct Part {
    float2 second,raw;
    float mass2,rawMass,baseCos;
};
inline Part transfer(float d,float eps,constant float2* nodes,uint order) {
    Part p={float2(0),float2(0),0,0,0};
    for(uint j=0;j<order;++j) {
        const float mag=11.0f*nodes[j].x,w2=mag*mag;
        const float weight=11.0f*nodes[j].y*invSqrt2Pi*sqrt(d)*exp(-.5f*d*w2);
        const float cw=cos(mag),sw=sin(mag);
        p.baseCos+=2.0f*weight*cw;
        for(uint i=0;i<2u;++i) {
            const float sg=i==0u?1.0f:-1.0f,w=sg*mag,x=eps*w;
            const float one=1.0f+x,s=sqrt(one);
            const float firstLog=w*(-2.5f+.5f*d*w2);
            const float secondLog=w2*(-2.5f*log2rem(x)-.5f*d*w2/one);
            const float secondPhase=w*w2*(s+3.0f)/(4.0f*(1.0f+s)*(1.0f+s)*(1.0f+s));
            const float2 y=float2(firstLog+eps*secondLog,-w2/((1.0f+s)*(1.0f+s)));
            const float2 e0=float2(cw,sg*sw);
            p.second+=weight*mul(e0,float2(secondLog,secondPhase)+exp2norm(y,eps));
            p.mass2+=weight*(secondLog+exp2real(y.x,eps));
            const float f=eps*y.y;
            const float scalar=weight*exp(eps*y.x);
            p.raw+=scalar*mul(e0,float2(cos(f),sin(f)));
            p.rawMass+=scalar;
        }
    }
    return p;
}
kernel void wells_scan(constant float2* nodes [[buffer(0)]],device Row* out [[buffer(1)]],
                       constant uint4& pars [[buffer(2)]],uint id [[thread_position_in_grid]],
                       uint lane [[thread_index_in_simdgroup]]) {
    const uint index=pars.z+id,g=id/32u;
    const uint band=(g/5u)%4u,mode=g%5u;
    const uint seed=hash32(index^0xa28ac144u^pars.w);
    const uint widths[4]={5u,8u,8u,5u},starts[4]={4u,9u,17u,25u};
    const uint exponent=starts[band]+hash32(seed+167169u)%widths[band];
    const uint base=1u<<exponent,n=base+(hash32(seed+687u)&(base-1u));
    const float turns=mode<4u?.25f*float(mode):unit(seed+76319u);
    const float eps=1.0f/(piValue*(float(n)+turns));
    const float c=mode==0u?1.0f:mode==2u?-1.0f:mode<4u?0.0f:cos(2.0f*piValue*turns);
    const float s=mode==1u?1.0f:mode==3u?-1.0f:mode<4u?0.0f:sin(2.0f*piValue*turns);
    const float angle=2.0f/eps,rawS=sin(angle),rawC=cos(angle);
    const Part x=transfer(1.0f,eps,nodes,pars.y);
    const Part y=transfer(2.0f,eps,nodes,pars.y);
    const float B=.25f,B2=-.25f*exp(-.25f);
    const float2 z=B*x.second+B2*y.second;
    const float M=.5f*(x.mass2+y.mass2);
    const float first=lead*s,second=M+c*z.x-s*z.y;
    const float stable=first+eps*second;
    const float2 raw=B*x.raw+B2*y.raw;
    const float rawScaled=(.5f*(x.rawMass+y.rawMass)-1.0f+rawC*raw.x-rawS*raw.y)/eps;
    const float quadFloor=B*x.baseCos+B2*y.baseCos;
    const bool good=isfinite(second)&&isfinite(rawScaled);
    for(uint metric=0;metric<3u;++metric) {
        float score=abs(rawS-s);
        if(metric==1u) score=abs(rawScaled-stable);
        if(metric==2u) score=abs(second-(.75f+c*secondCos));
        const float val=good?score:-INFINITY;
        const float mx=simd_max(val);
        const uint winner=simd_min(val==mx?lane:32u);
        const uint bad=simd_sum(uint(!good));
        if(lane==winner) out[g*3u+metric]=Row{eps,turns,s,rawS,first,second,stable,
             rawScaled,z.x,z.y,M,quadFloor,mx,n,mode,band,32u,index,bad,pars.w};
    }
}
