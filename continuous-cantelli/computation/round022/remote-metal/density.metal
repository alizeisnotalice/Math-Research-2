#include <metal_stdlib>
using namespace metal;

struct Core {
    float x;
    float c;
    float d;
    float radius;
    float variance;
};

struct Settings {
    uint shiftCount;
    uint epsilonCount;
    uint panelCount;
    float epsilonMin;
    float epsilonMax;
};

// Returns the local density divided by the target density, per unit u.
float relativeIntegrand(float u, float x, float shift, float variance, float targetLogScale, float targetHalfExponent) {
    float displacement = shift - u;
    float input = x + u;
    float exponent = targetLogScale + targetHalfExponent - 0.5f * (input * input + displacement * displacement / variance) - 0.5f * log(variance);
    return exp(exponent);
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
    float targetVariance = core.variance + epsilon;
    float target = core.x + shift;
    float targetHalfExponent = 0.5f * target * target / targetVariance;
    float targetLogScale = 0.5f * log(targetVariance / (2.0f * M_PI_F));
    float split = max(shift, 0.0f);
    float total = 0.0f;
    for (uint panel = 0; panel < settings.panelCount; ++panel) {
        // Left u = -r*s² resolves the linear cusp at zero.
        // On the right, quadratic distance from y resolves possible narrow peaks.
        float s = (float(panel) + 0.5f) / float(settings.panelCount);
        float s2 = s * s;
        float left = core.radius * s2;
        float lowRight = split * (1.0f - s2);
        float highRight = split + (core.radius - split) * s2;
        float leftTerm = relativeIntegrand(-left, core.x, shift, core.c * left + epsilon, targetLogScale, targetHalfExponent);
        float lowTerm = relativeIntegrand(lowRight, core.x, shift, core.d * lowRight * lowRight * lowRight + epsilon, targetLogScale, targetHalfExponent);
        float highTerm = relativeIntegrand(highRight, core.x, shift, core.d * highRight * highRight * highRight + epsilon, targetLogScale, targetHalfExponent);
        total += 2.0f * s * (core.radius * leftTerm + split * lowTerm + (core.radius - split) * highTerm);
    }
    ratios[index] = total / float(settings.panelCount);
}
