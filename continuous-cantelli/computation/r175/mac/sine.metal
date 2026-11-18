#include <metal_stdlib>
using namespace metal;
struct Settings {uint mode,seed,steps,unused;};
struct Candidate {float4 shape,origin,endpoint,sweep;};
uint mixbits(uint x) {
 x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16);
}
float unit(uint x) {return float(mixbits(x)>>8)*0x1p-24f;}
float2 sides(float u,float lambda) {
 if(lambda==0.f)return float2(u,1.f-u);
 return sin(lambda*float2(u,1.f-u));
}
kernel void positive_sine(constant Settings& cfg [[buffer(0)]],device Candidate *out [[buffer(1)]],
                          uint idx [[thread_position_in_grid]]) {
 const float X[8]={.0198550717512319f,.101666761293187f,.237233795041836f,.408282678752175f,
                   .591717321247825f,.762766204958164f,.898333238706813f,.980144928248768f};
 const float W[8]={.0506142681451881f,.111190517226687f,.156853322938944f,.181341891689181f,
                   .181341891689181f,.156853322938944f,.111190517226687f,.0506142681451881f};
 uint key=mixbits(idx+cfg.seed)+0xf123bac5u;
 float w1=exp2(-10.f+7.2f*unit(key)),w2=exp2(-10.f+7.2f*unit(key+1));
 float h1=.17f+1.08f*unit(key+2),h2=.17f+1.08f*unit(key+3);
 float c=3.5f*(unit(key+4)-.5f),L=exp2(-2.f+3.f*unit(key+5)),s=1.f;
 if(cfg.mode==0u){w2=w1;h2=h1;c=0.f;}
 if(cfg.mode==2u)s=.65f+1.05f*unit(key+6);
 if(cfg.mode==3u) {
    w1=.006f+.009f*unit(key);w2=w1;h1=.58f+.36f*unit(key+2);h2=h1;c=0.f;
    if(idx==0u){w1=.01f;w2=.01f;h1=.75f;h2=.75f;L=1.f;}
 }
 float3 nodes[48];
 float baseMass=0.f,sourceMass=0.f;
 float denom[8];
 for(uint i=0;i<8;i++){
    float y=c+L*(X[i]-.5f);
    denom[i]=W[i]*rsqrt(s)*exp(-.5f*y*y/s);
    baseMass+=denom[i];
 }
 uint cursor=0;
 for(uint side=0;side<2;side++){
    float w=side==0u?w1:w2;
    float h=side==0u?h1:h2;
    for(uint piece=0;piece<3;piece++){
       float begin=piece==0u?0.f:(piece==1u?.05f:.95f);
       float width=piece==1u?.90f:.05f;
       for(uint i=0;i<8;i++){
           float t=begin+width*X[i],u=side==0u?w*t:1.f-w*t;
           float x=c+L*(u-.5f);
           float mass=W[i]*width*w*exp(-.5f*x*x);
           float f=piece==0u?X[i]:(piece==1u?1.f:1.f-X[i]);
           nodes[cursor++]=float3(u,mass,h*f);
           sourceMass+=mass;
       }
    }
 }
 float2 linear=0,first=0,above=-1,pivot=0;
 for(uint j=0;j<cfg.steps;j++){
    float lambda=M_PI_F*float(j)/float(cfg.steps-1u),k2=lambda*lambda;
    float2 d=0,n=0;
    for(uint z=0;z<8;z++)d+=denom[z]*sides(X[z],lambda);
    for(uint z=0;z<48;z++){
       float3 r=nodes[z];
       n+=r.y*exp(k2*r.z)*sides(r.x,lambda);
    }
    float2 ratio=n/d;
    if(j==0u)linear=ratio;
    else if(j==cfg.steps-1u)first=ratio;
    else {
       if(ratio.x>above.x){above.x=ratio.x;pivot.x=lambda;}
       if(ratio.y>above.y){above.y=ratio.y;pivot.y=lambda;}
    }
 }
 out[idx]={float4(w1,w2,h1,h2),float4(c,L,s,0),float4(sourceMass/baseMass,linear,max(first.x,first.y)),
           float4(above,pivot)};
}
