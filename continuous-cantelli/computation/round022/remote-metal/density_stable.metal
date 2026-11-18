#include <metal_stdlib>
using namespace metal;

struct Core {
    float x;
    float c;
    float d;
    float radius;
    float excessVariance;
};

struct Settings {
    uint shiftCount;
    uint epsilonCount;
    uint panelCount;
    float epsilonMin;
    float epsilonMax;
};

float relativeIntegrand(float u, float x, float shift, float variance,
                        float targetExcess, float targetVariance, float targetLogScale) {
    float displacement = shift - u;
    float input = x + u;
    // Expand the difference of Gaussian squares before doing float32 arithmetic.
    float halfExponent = 0.5f * ((displacement * (2.0f * x + shift + u) - targetExcess * input * input)
                               / targetVariance - displacement * displacement / variance);
    return exp(targetLogScale + halfExponent - 0.5f * log(variance));
}

kernel void scan_density(
    device const Core* cores [[buffer(0)]],
    constant Settings& settings [[buffer(1)]],
    device float* ratios [[buffer(2)]],
    uint index [[thread_position_in_grid]]) {
    uint pointsPerCore = settings.shiftCount * settings.epsilonCount;
    uint coreIndex = index / pointsPerCore;
    uint shiftIndex = (index % pointsPerCore) / settings.epsilonCount;
    uint epsilonIndex = index % settings.epsilonCount;
    Core core = cores[coreIndex];
    float shift = core.radius * (2.0f * float(shiftIndex) / float(settings.shiftCount - 1) - 1.0f);
    float epsilonWeight = float(epsilonIndex) / float(settings.epsilonCount - 1);
    float epsilon = exp(mix(log(settings.epsilonMin), log(settings.epsilonMax), epsilonWeight));
    float targetExcess = core.excessVariance + epsilon;
    float targetVariance = 1.0f + targetExcess;
    float targetLogScale = 0.5f * log(targetVariance / (2.0f * M_PI_F));
    float split = max(shift, 0.0f);
    float total = 0.0f;
    for (uint panel = 0; panel < settings.panelCount; ++panel) {
        // Quadratic substitutions resolve the linear cusp and narrow right peak.
        float s = (float(panel) + 0.5f) / float(settings.panelCount);
        float s2 = s * s;
        float left = core.radius * s2;
        float lowRight = split * (1.0f - s2);
        float highRight = split + (core.radius - split) * s2;
        float leftTerm = relativeIntegrand(-left, core.x, shift, core.c * left + epsilon, targetExcess, targetVariance, targetLogScale);
        float lowTerm = relativeIntegrand(lowRight, core.x, shift, core.d * lowRight * lowRight * lowRight + epsilon, targetExcess, targetVariance, targetLogScale);
        float highTerm = relativeIntegrand(highRight, core.x, shift, core.d * highRight * highRight * highRight + epsilon, targetExcess, targetVariance, targetLogScale);
        total += 2.0f * s * (core.radius * leftTerm + split * lowTerm + (core.radius - split) * highTerm);
    }
    ratios[index] = total / float(settings.panelCount);
}
