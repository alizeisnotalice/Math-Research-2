#include <metal_stdlib>
using namespace metal;
struct Probe {float logK,logE;uint nearOne,pad;};
float logSum(float x,float k) {
 float term=2.f*x,sum=term,shift=0.f;
 for(uint n=1;n<1600u;n++){
    float f=float(n);
    term*=((f-k)*x)/((f+.5f)*(f+1.f));
    sum+=term;
    if(sum>0x1p60f) {sum*=0x1p-60f;term*=0x1p-60f;shift+=41.58883083359672f;}
    if(f>x+8.f && term<2.e-8f*sum)break;
 }
 return log(sum)+shift;
}
float higher(float x,float k) {
 float term=(2.f/3.f)*x*x,sum=term;
 for(uint n=2;n<48u;n++){
    float f=float(n);
    term*=((f-k)*x)/((f+.5f)*(f+1.f));
    sum+=term;
    if(term<1.e-9f*sum)break;
 }
 return sum;
}
kernel void marchaud(device const Probe* probes [[buffer(0)]],device float4* out [[buffer(1)]],
                      uint idx [[thread_position_in_grid]]) {
 Probe p=probes[idx];
 if (p.nearOne!=0u) {
     float e=exp(p.logE),k=1.f-e;
     float lo=0.f,hi=.5f/k;
     for(uint j=0;j<25u;j++){
         float r=.5f*(lo+hi),x=fma(e,r,.5f);
         float f=fma(2.f*k,r,k*higher(x,k))-1.f;
         if(f>0)hi=r;else lo=r;
     }
     float r=.5f*(lo+hi),x=fma(e,r,.5f),y=sqrt(2.f*x);
     out[idx]=float4(y,x,2.f*r/(y+1.f),fma(2.f*k,r,k*higher(x,k))-1.f);
 } else {
     float k=exp(p.logK),ell=-p.logK;
     float lo=.5f,hi=max(2.f,ell+.5f*log(max(1.f,ell))+5.f);
     for(uint j=0;j<3u && logSum(hi,k)+p.logK<=0.f;j++)hi*=2.f;
     for(uint j=0;j<25u;j++){
        float x=.5f*(lo+hi);
        if(logSum(x,k)+p.logK>0)hi=x;else lo=x;
     }
     float x=.5f*(lo+hi),y=sqrt(2.f*x);
     out[idx]=float4(y,x,(y-1.f)/(1.f-k),logSum(x,k)+p.logK);
 }
}
