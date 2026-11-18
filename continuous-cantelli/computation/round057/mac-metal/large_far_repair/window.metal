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
    // The small-shift branch is already covered by the finite scan.
    if (fabs(s)<=4.0f) {
        ratios[index]=as_type<float>(0x7f7fffffu);
        return;
    }
    float ratio=1.0f+s/p.a;
    float q=ratio*ratio;
    if (q==0.0f) {
        ratios[index]=0.0f;
        return;
    }
    float gain=q*p.logGain;
    float logPositive=gain+log(-em1(-gain))-log(A);
    float logNegative=log(-em1(q*logLoss))-log(p.b);
    float observed=fabs(s);
    float dv=max(0.0f,observed-v);
    float db=observed-clamp(observed,p.r,p.r+p.length);
    float reference=max(logPositive-0.5f*dv*dv,logNegative-0.5f*db*db);
    float valley=0.0f,band=0.0f;
    for(uint i=0;i<settings.panels;++i) {
        float f=(float(i)+0.5f)/float(settings.panels);
        float yv=v*f;
        float yb=p.r+p.length*f;
        float sv0=s-yv,sv1=s+yv,sb0=s-yb,sb1=s+yb;
        valley+=exp(logPositive-reference-0.5f*sv0*sv0)+exp(logPositive-reference-0.5f*sv1*sv1);
        band+=exp(logNegative-reference-0.5f*sb0*sb0)+exp(logNegative-reference-0.5f*sb1*sb1);
    }
    float norm=1.0f/(sqrt(2.0f*M_PI_F)*float(settings.panels));
    float signedMass=norm*(v*valley/p.valleyMass-p.length*band/p.bandMass);
    // Retain a checkable sign even if the physical residual underflows.
    ratios[index]=signedMass*exp(max(reference,-70.0f))*(p.a*p.a*p.moment/(derivative*derivative));
}
