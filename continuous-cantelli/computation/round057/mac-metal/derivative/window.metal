#include <metal_stdlib>
using namespace metal;

struct Core { float a,w,r,length,b,moment,bandMass,valleyMass,logGain; };
struct Settings { uint shifts,panels; };

float em1(float x) {
    if (fabs(x)>0.10f) return exp(x)-1.0f;
    return x*(1.0f+x*(0.5f+x*(0.16666666667f+x*(0.04166666667f+x*(0.00833333333f+x*0.00138888889f)))));
}
float shapeDelta(float y,float s,float shiftDelta,float shiftFactor) {
    float sy=fabs(s*y);
    float e=em1(sy);
    return exp(-0.5f*y*y)*(shiftDelta+shiftFactor*(0.5f*e*e*exp(-sy)));
}

kernel void window_scan(device const Core* cores [[buffer(0)]],
                      constant Settings& settings [[buffer(1)]],
                      device float* ratios [[buffer(2)]],
                      uint index [[thread_position_in_grid]]) {
    Core p=cores[index/settings.shifts];
    float nu=-2.0f+10.0f*float(index%settings.shifts)/float(settings.shifts-1);
    float A=p.b*p.bandMass/p.valleyMass;
    float logLoss=log(1.0f-p.b);
    float derivative=(1.0f+1.0f/A)*p.logGain+(1.0f-p.b)/p.b*logLoss;
    float s=nu*derivative/(p.a*p.moment);
    float v=p.w/p.a;
    float shiftDelta=em1(-0.5f*s*s);
    float shiftFactor=exp(-0.5f*s*s);
    float aValley=0.0f,aBand=0.0f;
    for(uint i=0;i<settings.panels;++i) {
        float f=(float(i)+0.5f)/float(settings.panels);
        aValley+=shapeDelta(v*f,s,shiftDelta,shiftFactor);
        aBand+=shapeDelta(p.r+p.length*f,s,shiftDelta,shiftFactor);
    }
    float norm=2.0f/(sqrt(2.0f*M_PI_F)*float(settings.panels));
    float valleyChange=norm*v*aValley/p.valleyMass;
    float bandChange=norm*p.length*aBand/p.bandMass;
    float qa=s/p.a;
    float dq=qa*(2.0f+qa);
    float aChange=(1.0f+1.0f/A)*em1(dq*p.logGain);
    float bChange=-(1.0f-p.b)/p.b*em1(dq*logLoss);
    float res=(aChange+valleyChange+aChange*valleyChange)-(bChange+bandChange+bChange*bandChange);
    ratios[index]=res*(p.a*p.a*p.moment/(derivative*derivative));
}
