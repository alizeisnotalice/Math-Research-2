#include <metal_stdlib>
using namespace metal;
struct Node {float u,re,im,pad;};
struct Config {float step,width,offset,scaleM;uint m,mid,order,pad;};
kernel void cantor_window(device const Node *nodes [[buffer(0)]],
                          device const float2 *scales [[buffer(1)]],
                          device const float2 *rotations [[buffer(2)]],
                          constant Config &cfg [[buffer(3)]],
                          device float4* out [[buffer(4)]],
                          uint idx [[thread_position_in_grid]]) {
 float w=(float(idx)-float(cfg.mid))*cfg.step;
 float2 pref[13];
 for(uint j=0;j<cfg.m;j++){
    float2 ab=scales[j];
    float angle=w*ab.x;
    float delta=fma(w,ab.x,-angle)+w*ab.y;
    float sa=sin(angle),ca=cos(angle),sd=sin(delta),cd=cos(delta);
    pref[j]=float2(ca*cd-sa*sd,sa*cd+ca*sd);
 }
 float2 base=float2(0),resonant=float2(0);
 for(uint q=0;q<cfg.order;q++){
    float factor=1;
    for(uint j=0;j<cfg.m;j++){
       float2 ru=rotations[q*13u+j];
       factor*=pref[j].x*ru.x-pref[j].y*ru.y;
    }
    Node node=nodes[q];
    float t=fma(node.u,cfg.scaleM,w*cfg.scaleM);
    float2 x=float2(t,t+2.f*M_PI_F),tail=float2(1);
    for(uint j=0;j<12u;j++){
      x*=1.f/3.f;
      tail*=cos(x);
    }
    float2 z=factor*float2(node.re,node.im);
    base+=z*tail.x;
    resonant+=z*tail.y;
 }
 float r=dot(resonant,resonant),b=dot(base,base);
 float y=w/cfg.width,d=(w-cfg.offset)/cfg.width;
 out[idx]=float4(r*exp(-.5f*d*d),b*exp(-.5f*y*y),r,b);
}
