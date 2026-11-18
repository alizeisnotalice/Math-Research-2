#include <metal_stdlib>
using namespace metal;

struct Candidate {
    float m;
    float r;
    float epsilon;
    float delta;
};

// On the first panel set u = b*t^power to resolve the low-variance contact.
inline float panel(Candidate c, constant float2* nodes, uint n, float a, float b,
                   uint power, bool right, float exponentOffset, float prefactor) {
    if (b <= a) return 0.0f;
    float result = 0.0f;
    float width = b - a;
    for (uint i = 0; i < n; ++i) {
        float t = nodes[i].x;
        float t2 = t * t;
        float u = a + width * (power == 1 ? t : power == 2 ? t2 : t2 * t);
        float jacobian = width * (power == 1 ? 1.0f : power == 2 ? 2.0f * t : 3.0f * t2);
        float x = right ? 1.0f + u : 1.0f - u;
        float z = right ? c.delta - u : c.delta + u;
        float variance = c.epsilon + (right ? u * u * u : u);
        float exponent = exponentOffset - 0.5f * (x*x + z*z / variance);
        result += (nodes[i].y * jacobian * prefactor) * rsqrt(variance) * exp(exponent);
    }
    return result;
}

kernel void core_ratio(device const Candidate* grid [[buffer(0)]],
                       constant float2* nodes [[buffer(1)]],
                       device float* output [[buffer(2)]],
                       constant uint2& counts [[buffer(3)]],
                       uint index [[thread_position_in_grid]]) {
    if (index >= counts.x) return;
    Candidate c = grid[index];
    float targetVariance = c.m + c.epsilon;
    float y = 1.0f + c.delta;
    float exponentOffset = 0.5f * y*y / targetVariance;
    float prefactor = sqrt(targetVariance * 0.15915494309189535f);
    float left = panel(c, nodes, counts.y, 0.0f, c.r, 2, false, exponentOffset, prefactor);
    float right;
    if (c.delta > 0.0f) {
        float s = sqrt(c.epsilon + c.delta*c.delta*c.delta);
        float a = clamp(c.delta - 8.0f * s, 0.0f, c.r);
        float b = clamp(c.delta + 8.0f * s, 0.0f, c.r);
        right = panel(c, nodes, counts.y, 0.0f, a, 3, true, exponentOffset, prefactor)
              + panel(c, nodes, counts.y, a, b, a == 0.0f ? 3 : 1, true, exponentOffset, prefactor)
              + panel(c, nodes, counts.y, b, c.r, b == 0.0f ? 3 : 1, true, exponentOffset, prefactor);
    } else {
        float b = min(8.0f * sqrt(c.epsilon), c.r);
        right = panel(c, nodes, counts.y, 0.0f, b, 3, true, exponentOffset, prefactor)
              + panel(c, nodes, counts.y, b, c.r, 3, true, exponentOffset, prefactor);
    }
    output[index] = left + right;
}
