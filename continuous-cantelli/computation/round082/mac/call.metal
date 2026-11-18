#include <metal_stdlib>
using namespace metal;
struct Core {float k,l,coefficient,lambda,center,mean,std;};
struct Settings {uint panels,strikes;};
kernel void scan_calls(device const Core* cores [[buffer(0)]],
                      constant Settings& cfg [[buffer(1)]],
                      device float* result [[buffer(2)]],
                      uint idx [[thread_position_in_grid]]) {
    Core p=cores[idx/cfg.strikes];
    float y=7.0f*float(idx%cfg.strikes)/float(cfg.strikes-1);
    float target=sqrt(1.0f+p.mean);
    float strike=y*target;
    const float4 u=float4(0.0694318442f,0.3300094782f,0.6699905218f,0.9305681558f);
    const float4 w=float4(0.1739274226f,0.3260725774f,0.3260725774f,0.1739274226f);
    float total=0.0f;
    for(uint j=0;j<cfg.panels;++j) {
        float x=9.0f*(float(j)+0.5f)/float(cfg.panels);
        float change=p.lambda*(cos(p.k*x)+p.coefficient*cos(p.l*x)-p.center);
        float4 noise=p.mean+u*change;
        float4 inv=1.0f/noise;
        float minus=strike-x,plus=strike+x;
        float4 density=(exp((-0.5f*minus*minus)*inv)+exp((-0.5f*plus*plus)*inv))*sqrt(inv);
        total+=exp(-0.5f*x*x)*change*dot(w,density);
    }
    result[idx]=total*9.0f/(4.0f*M_PI_F*float(cfg.panels)*target);
}
