#include <metal_stdlib>
using namespace metal;
struct Settings {uint seed,mode,count,unused;};
struct Candidate {float4 point,stable,raw,coeff;};
uint mixbits(uint x) {
 x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);
}
float unit(uint x) {return float(mixbits(x)>>8)*0x1p-24f;}
float tSecond(float v) {
 return v-v*v+(3.f/8.f)*(3.f*v*v-v*v*v)-(3.f*v*v-6.f*v*v*v+v*v*v*v)/32.f;
}
kernel void two_wells(constant Settings& cfg [[buffer(0)]],device Candidate *out [[buffer(1)]],
                      uint idx [[thread_position_in_grid]]) {
 const float b1=.5f,b3=-.3582656552868946f;
 const float c1=.6065306597126334f,c3=.8464817248906141f;
 const uint seed=mixbits(idx+cfg.seed)+0x9e3779b9u;
 const uint panels=16u<<(cfg.mode%6u), phaseClass=cfg.mode/6u;
 uint k=uint(round(exp2(1.f+15.f*unit(seed))));
 float cls=float(phaseClass),phi=cls*(M_PI_F*.5f);
 float co=0,si=0;
 if (phaseClass==0u)co=1;
 else if (phaseClass==1u)si=1;
 else if (phaseClass==2u)co=-1;
 else if (phaseClass==3u)si=-1;
 else {phi=2.f*M_PI_F*unit(seed+1);co=cos(phi);si=sin(phi);}
 float R=M_PI_F*float(k)+.5f*phi;
 float invR=1.f/R, rr=invR*invR;
 float2 leading=float2(1+b1*c1*co,1+b3*c3*co);
 float2 r1=float2(b1*c1*si*1.5f,b3*c3*si*(4.f/9.f));
 float2 r2=float2(1.f+b1*c1*co*tSecond(1.f),1.f/3.f+b3*c3*co*tSecond(1.f/3.f));
 float step=24.f/float(panels);
 float2 raw=0,res=0;
 for (uint i=0;i<panels/2u;i++) {
    float y=(float(i)+.5f)*step, h=.5f*y*invR;
    float v=sqrt(1.f+h*h);
    float evenDelta=4.f*h*h, evenW=1.f+evenDelta;
    float oddW=-h*(3.f+4.f*h*h)/v;
    float eta=.5f*y*y*invR/(v+1.f);
    float sy=sin(y),cy=cos(y), se=sin(eta),dce=-2.f*sin(.5f*eta)*sin(.5f*eta);
    float shiftCo=co*dce-si*se;
    float shiftedSi=si*(1.f+dce)+co*se;
    float correctionOsc=evenDelta*co*cy+evenW*cy*shiftCo-oddW*sy*shiftedSi;
    float2 g=exp(float2(-.5f*y*y,-1.5f*y*y))*float2(.3989422804014327f,.690988298942671f);
    float2 bs=float2(b1,b3);
    res+=g*(evenDelta+bs*correctionOsc);
    float t1=co*cy+shiftCo*cy, t2=shiftedSi*sy;
    raw+=g*(evenW*(1.f+bs*t1)-bs*oddW*t2);
 }
 raw*=2.f*step;res*=2.f*step;
 float pair=res.x+res.y;
 float residual2=pair-(r1.x+r1.y)*invR-(r2.x+r2.y)*rr;
 out[idx]={float4(float(k),phi,R,float(panels)),float4(res,pair,residual2),
         float4(raw,raw.x+raw.y-2.f,max(abs(raw.x-leading.x-res.x),abs(raw.y-leading.y-res.y))),
         float4(r1.x+r1.y,r2.x+r2.y,leading.x,leading.y)};
}
