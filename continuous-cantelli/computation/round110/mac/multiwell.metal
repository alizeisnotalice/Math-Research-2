#include <metal_stdlib>
using namespace metal;

struct Settings {uint seed, panels, mode, count;};
struct Candidate {float4 roots; float4 diagnostics; float4 margins;};

uint mixbits(uint x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15;
    x *= 0x846ca68bu; return x ^ (x >> 16);
}
float uniform01(uint x) {return float(mixbits(x) >> 8) * 0x1p-24f;}
float poly(float s, float4 roots) {
    return s*(s-roots.x)*(s-roots.y)*(s-roots.z)*(s-roots.w);
}
float minimumPoint(float lo, float hi, float4 roots) {
    for (uint k=0; k<17; ++k) {
        float m=(lo+hi)*0.5f;
        float4 d=1.0f/(m-roots);
        float slope=1.0f/m+d.x+d.y+d.z+d.w;
        if (slope>0.0f) lo=m; else hi=m;
    }
    return (lo+hi)*0.5f;
}
float remainder(float v,float w) {
    if (fabs(w)<0.04f)
        return exp(v)*w*w*(0.5f+w*(1.0f/6.0f+w*(1.0f/24.0f+w/120.0f)));
    return exp(v+w)-exp(v)*(1.0f+w);
}
kernel void multiwell(constant Settings& cfg [[buffer(0)]],
                      device Candidate *result [[buffer(1)]],
                      uint idx [[thread_position_in_grid]]) {
    uint seed=mixbits(idx+cfg.seed)+0x9e3779b9u;
    float u=exp2(-3.0f+7.0f*uniform01(seed));
    float r1=.02f+.49f*uniform01(seed+1);
    float r2=.02f+.73f*uniform01(seed+2);
    float r3=.02f+.96f*uniform01(seed+3);
    if (cfg.mode==1) {
        r1=.02f+.24f*uniform01(seed+1);
        r2=.27f+.32f*uniform01(seed+2);
        r3=.60f+.37f*uniform01(seed+3);
    }
    if (cfg.mode==2) {
        u=exp2(-2.0f+5.0f*uniform01(seed));
        r1=.05f+.30f*uniform01(seed+1);
        r2=.22f+.50f*uniform01(seed+2);
        r3=.43f+.54f*uniform01(seed+3);
    }
    float low=min(r1,min(r2,r3)), high=max(r1,max(r2,r3));
    float middle=r1+r2+r3-low-high;
    r1=low; r2=middle; r3=high;
    Candidate out={float4(r1,r2,r3,u),float4(-1),float4(0)};
    result[idx]=out;
    float A=r1+r2+r3, B=r1*r2+r1*r3+r2*r3, C=r1*r2*r3;
    float g[16];
    for (uint j=0; j<16; ++j) g[j]=rsqrt(1.0f+u*float(j));
    float den=g[4]-A*g[3]+B*g[2]-C*g[1];
    if (fabs(den)<2e-5f) return;
    float r4=(g[5]-A*g[4]+B*g[3]-C*g[2])/den;
    if (!(r4>0.025f && r4<0.975f) || min(min(r2-r1,r3-r2),min(abs(r4-r1),min(abs(r4-r2),abs(r4-r3))))<0.023f) return;
    float p[6]={0,C*r4,-C-B*r4,B+A*r4,-A-r4,1};
    float h2=0, norm=0;
    for (uint i=1; i<=5; ++i) {
        h2 += p[i]*(-float(i)*u)*g[i]*g[i]*g[i];
        for (uint j=1; j<=5; ++j) norm += p[i]*p[j]*g[i+j];
    }
    if (!(h2>0 && norm>1e-6f)) return;
    float eps=2.0f*h2/norm;
    if (!(eps>0.025f && eps<20000.0f)) return;
    float4 roots=float4(min(r1,r4),max(r1,min(r2,r4)),max(r2,min(r3,r4)),max(r3,r4));
    float z0=minimumPoint(roots.x,roots.y,roots);
    float z1=minimumPoint(roots.z,roots.w,roots);
    float p0=poly(z0,roots),p1=poly(z1,roots);
    float fmin=1.0f+eps*min(p0,p1);
    if (!(fmin>0.03f && fmin<0.99999f)) return;
    float q=fmin-u/(1.0f+u);
    float xstar=sqrt(-2.0f*log(p0<=p1?z0:z1)/u);
    float score=0, jensen=1.0e20f, minGap=1.0e20f, maxGap=-1.0e20f;
    float mask=0, old=0;
    const float lambdas[12]={.25f,.5f,.75f,1.f,1.5f,2.f,2.5f,3.f,4.f,6.f,8.f,12.f};
    for (uint j=0; j<12; ++j) {
        float l=lambdas[j], t=0.5f*l*l;
        float shift=0;
        for (uint k=1; k<=5; ++k)
            shift += p[k]*g[k]*exp(-t*(float(k)*u/(1.0f+float(k)*u)-u/(1.0f+u)));
        jensen=min(jensen,shift);
        float nonlinear=0;
        for (uint k=0; k<cfg.panels; ++k) {
            float y=28.0f*(float(k)+0.5f)/float(cfg.panels)-14.0f;
            float x=l+y;
            float s=exp(-0.5f*u*x*x);
            float v=-0.5f*y*y+t*u/(1.0f+u);
            float w=-t*eps*poly(s,roots);
            nonlinear+=remainder(v,w);
        }
        float normalized=-shift+0.3989422804f*28.0f*nonlinear/(float(cfg.panels)*t*eps);
        minGap=min(minGap,normalized);
        maxGap=max(maxGap,normalized);
        score=max(score,abs(normalized));
        if (j>0 && old*normalized<0.0f) mask+=1;
        old=normalized;
    }
    out.diagnostics=float4(eps,fmin,q,xstar);
    out.margins=float4(score,jensen,minGap,mask);
    result[idx]=out;
}
