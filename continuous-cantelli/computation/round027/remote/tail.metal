#include <metal_stdlib>
using namespace metal;

struct Profile { float lambda, h, m, downWidth, plateauWidth, logBase; };
struct Settings { uint offsets, panels; };

float smoothEdge(float s) {
    float left = exp(-1.0f/s);
    float right = exp(-1.0f/(1.0f-s));
    return left / (left+right);
}

float scaledIntegrand(float t, float q, float offset, float invScale2, float reference, float penalty) {
    float z = offset + t;
    return exp(-t - 0.5f*z*z*invScale2 - penalty*q - reference) * (1.0f-exp(-penalty*(1.0f-q)));
}

kernel void scan_tail(device const Profile* profiles [[buffer(0)]],
                      constant Settings& settings [[buffer(1)]],
                      device float2* out [[buffer(2)]],
                      uint index [[thread_position_in_grid]]) {
    Profile p = profiles[index / settings.offsets];
    float kappa = -4.0f + 16.0f * float(index % settings.offsets) / float(settings.offsets-1);
    float logLambda = log(p.lambda);
    float offset = kappa*logLambda;
    float l2 = p.lambda*p.lambda;
    float invScale2 = 1.0f/(l2*(1.0f-p.m));
    float penalty = l2*(p.h-p.m)/2.0f;
    float reference = p.downWidth/2.0f;
    float total = 0.0f;
    for (uint i=0; i<settings.panels; ++i) {
        float s = (float(i)+0.5f)/float(settings.panels);
        float ramp = smoothEdge(s);
        total += p.downWidth*scaledIntegrand(-p.downWidth*s,ramp,offset,invScale2,reference,penalty);
        total += p.plateauWidth*scaledIntegrand(p.plateauWidth*s,0.0f,offset,invScale2,reference,penalty);
        total += p.downWidth*scaledIntegrand(p.plateauWidth+p.downWidth*s,ramp,offset,invScale2,reference,penalty);
    }
    float scale = p.lambda*sqrt(1.0f-p.m)*sqrt(2.0f*M_PI_F);
    float logCorrection = log(total/float(settings.panels))-offset+reference-log(scale);
    float m = max(logCorrection,p.logBase);
    out[index] = float2(m+log(1.0f+exp(min(logCorrection,p.logBase)-m)),logCorrection);
}
