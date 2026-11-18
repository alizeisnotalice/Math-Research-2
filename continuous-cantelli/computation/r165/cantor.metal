#include <metal_stdlib>
using namespace metal;
struct Level {float length,shift,tail,pad;};
struct Case {float variance,c;uint geometry,pad;};
struct Settings {uint points,begin,uniformPoints,pad;};
float address_center(device const Level *g,uint level,uint addr) {
 float center=0;
 for(uint k=1;k<=level;k++)
    center+=((addr>>(level-k))&1u)?g[k].shift:-g[k].shift;
 return center;
}
uint lower_index(device const Level *g,uint level,float threshold) {
 uint index=0;float center=0;
 for(uint k=1;k<=level;k++){
    bool right=threshold>center;
    index=(index<<1u)|uint(right);
    center+=right?g[k].shift:-g[k].shift;
 }
 return index+uint(center<threshold);
}
kernel void cantor_scan(device const Level *geom [[buffer(0)]],
                        device const Case *cases [[buffer(1)]],
                        constant Settings& cfg [[buffer(2)]],
                        device float4* out [[buffer(3)]],
                        uint idx [[thread_position_in_grid]]) {
 uint caseIndex=cfg.begin+idx/cfg.points, site=idx%cfg.points;
 Case test=cases[caseIndex];
 device const Level *g=geom+test.geometry*22u;
 float sigma=sqrt(test.variance);
 float x;
 if(site<cfg.uniformPoints) x=.5f*float(site)/float(cfg.uniformPoints-1u);
 else {
    uint depth=1;
    while(depth<17u && g[depth].length>1.5f*sigma)depth++;
    uint n=1u<<(depth-1u),s=site-cfg.uniformPoints,budget=cfg.points-cfg.uniformPoints;
    uint addr=n+s%n;
    float off=0;
    if(n>=budget)addr=n+uint(floor(float(n)*((float(s)+.5f)/float(budget))));
    else {
       uint refinements=budget/n;
       off=((float(s/n)+.5f)/float(refinements)-.5f)*g[depth].length;
       if(s/n>=refinements)off=0;
    }
    x=address_center(g,depth,addr)+off;
 }
 uint depth=6u;
 while(depth<20u && g[depth].length>sigma/12.f)depth++;
 float width=g[depth].length, norm=width*g[depth].tail/(sigma*sqrt(2.f*M_PI_F));
 uint lo=lower_index(g,depth,x-8.f*sigma-.5f*width);
 uint hi=lower_index(g,depth,x+8.f*sigma+.5f*width);
 uint cap=1u<<depth;
 lo=min(lo,cap);hi=min(hi,cap);
 float delta=(width/sigma)*.3872983346207417f;
 float4 sum=0;
 for(uint a=lo;a<hi;a++) {
    float z=(address_center(g,depth,a)-x)/sigma;
    float3 u=float3(z-delta,z,z+delta);
    float3 w=exp(-.5f*u*u)*float3(5.f/18.f,4.f/9.f,5.f/18.f);
    sum.x+=w.x+w.y+w.z;
    sum.y+=dot(w,u);
    sum.z+=dot(w,u*u-1.f);
 }
 out[idx]=float4(norm*sum.xyz,x);
}
