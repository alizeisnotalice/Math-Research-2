#include <metal_stdlib>
using namespace metal;

struct Core {
    float center;
    float coefficient;
    float radius;
    float excessVariance;
};

struct Settings {
    uint shiftCount;
    uint epsilonCount;
    uint panels;
};

float relativeDensity(float u, float shift, float epsilon, Core core) {
    float u2 = u * u;
    float local = core.coefficient * u2 * u2 + epsilon;
    float excess = core.excessVariance + epsilon;
    float target = 1.0f + excess;
    float displacement = shift - u;
    float input = core.center + u;
    float exponent = 0.5f * ((displacement * (2.0f * core.center + shift + u) - excess * input * input) / target
                             - displacement * displacement / local);
    return exp(exponent) * sqrt(target / (2.0f * M_PI_F * local));
}

kernel void quartic_scan(device const Core* cores [[buffer(0)]],
                         constant Settings& settings [[buffer(1)]],
                         device float* ratios [[buffer(2)]],
                         uint index [[thread_position_in_grid]]) {
    uint pointsPerCore = settings.shiftCount * settings.epsilonCount;
    uint coreIndex = index / pointsPerCore;
    uint shiftIndex = (index % pointsPerCore) / settings.epsilonCount;
    uint epsilonIndex = index % settings.epsilonCount;
    Core core = cores[coreIndex];
    float shift = 1.25f * core.radius * (2.0f * float(shiftIndex) / float(settings.shiftCount - 1) - 1.0f);
    float epsilon = exp(mix(log(1e-7f), log(0.25f), float(epsilonIndex) / float(settings.epsilonCount - 1)));
    float split = clamp(shift, -core.radius, core.radius);
    float lowWidth = split + core.radius;
    float highWidth = core.radius - split;
    float total = 0.0f;
    // Quadratic distance from closest available input resolves the small-noise peak.
    for (uint panel = 0; panel < settings.panels; ++panel) {
        float s = (float(panel) + 0.5f) / float(settings.panels);
        float low = split - lowWidth * s * s;
        float high = split + highWidth * s * s;
        total += 2.0f * s * (lowWidth * relativeDensity(low, shift, epsilon, core) + highWidth * relativeDensity(high, shift, epsilon, core));
    }
    ratios[index] = total / float(settings.panels);
}
