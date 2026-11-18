#include <metal_stdlib>
using namespace metal;
struct Core { float a,w,r,length,b,moment,bandMass,valleyMass,logGain; };
struct Settings { uint shifts,panels; };

float e1(float x) {
    if(fabs(x)>0.10f) return exp(x)-1.0f;
    return x*(1.0f+x*(0.5f+x*(0.16666666667f+x*(0.04166666667f+x*0.00833333333f))));
}

kernel void window_scan(device const Core* cores [[buffer(0)]],
                        constant Settings& settings [[buffer(1)]],
                        device float* scores [[buffer(2)]],
                        uint index [[thread_position_in_grid]]) {
    Core p=cores[index/settings.shifts];
    float factor=0.5f+float(index%settings.shifts)/float(settings.shifts-1);
    float radius=p.r+p.length;
    float v=p.w/p.a;
    float s=factor*radius*p.a*p.a/p.logGain;
    float q=1.0f+s/p.a;
    q*=q;
    float gain=q*p.logGain;
    float loss=q*log(1.0f-p.b);
    float upperV=min(s*v,48.0f);
    float upperB=min(s*p.length,48.0f);
    float valley=0.0f,band=0.0f;
    for(uint i=0;i<settings.panels;++i) {
        float u=(float(i)+0.5f)/float(settings.panels);
        float zv=upperV*u,zb=upperB*u;
        float dv=zv/s,db=zb/s;
        float yv=v-dv,yb=radius-db;
        valley+=exp(-zv+dv*(v-0.5f*dv))*(1.0f+exp(-2.0f*s*yv));
        band+=exp(-zb+db*(radius-0.5f*db))*(1.0f+exp(-2.0f*s*yb));
    }
    float logCore=log(upperV*valley)-log(upperB*band);
    float logAmplitude=gain+log(-e1(-gain))-log(-e1(loss));
    float action=logAmplitude-s*(radius-v)+0.5f*(radius*radius-v*v)+logCore;
    scores[index]=action/(s*radius);
}
