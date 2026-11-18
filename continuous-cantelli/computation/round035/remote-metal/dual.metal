#include <metal_stdlib>
using namespace metal;

struct Core { float center, width, power, depth; };
struct Settings { uint offsets, halfWidths, panels; };

float gaussianTail(float z) {
    // Stable erfc approximation on the positive half-axis.
    float t = 1.0f/(1.0f+0.3535533905932738f*z);
    float p = 0.17087277f;
    p = -0.82215223f+t*p;
    p = 1.48851587f+t*p;
    p = -1.13520398f+t*p;
    p = 0.27886807f+t*p;
    p = -0.18628806f+t*p;
    p = 0.09678418f+t*p;
    p = 0.37409196f+t*p;
    p = 1.00002368f+t*p;
    return 0.5f*t*exp(-0.5f*z*z-1.26551223f+t*p);
}

float probability(float lower, float upper) {
    if (lower>=0) return gaussianTail(lower)-gaussianTail(upper);
    if (upper<=0) return gaussianTail(-upper)-gaussianTail(-lower);
    return 1.0f-gaussianTail(-lower)-gaussianTail(upper);
}

kernel void dual_scan(device const Core* cores [[buffer(0)]],
                      constant Settings& settings [[buffer(1)]],
                      device float* ratios [[buffer(2)]],
                      uint index [[thread_position_in_grid]]) {
    uint points = settings.offsets*settings.halfWidths;
    Core core = cores[index / points];
    uint offsetIndex = index % points / settings.halfWidths;
    uint halfIndex = index % settings.halfWidths;
    float offset = -4.0f+8.0f*float(offsetIndex)/float(settings.offsets-1);
    float halfWidth = exp(mix(log(0.03125f),log(4.0f),float(halfIndex)/float(settings.halfWidths-1)));
    float total=0.0f;
    for (uint i=0;i<settings.panels;++i) {
        float unit=-1.0f+2.0f*(float(i)+0.5f)/float(settings.panels);
        float y=unit*core.width;
        float gap=core.depth*pow(1.0f-unit*unit,core.power)/(core.center*core.center);
        float f=1.0f-gap;
        float x=core.center+y;
        float meanOffset=(y+core.center*gap)/f;
        float sigmaInverse=sqrt(f);
        float lo=(offset-halfWidth-meanOffset)*sigmaInverse;
        float hi=(offset+halfWidth-meanOffset)*sigmaInverse;
        float weight=exp(0.5f*x*x*gap/f)/sigmaInverse;
        total += weight*probability(lo,hi);
    }
    ratios[index]=total*core.width/(halfWidth*float(settings.panels));
}
