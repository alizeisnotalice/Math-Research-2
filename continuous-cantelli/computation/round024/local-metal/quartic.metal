#include <metal_stdlib>
using namespace metal;

struct Point {
    float m, a, c, r, epsilon, delta, margin;
};
struct Row {
    float ratio, m, a, c, r, epsilon, delta, margin;
    uint tested, above, nonfinite, index;
};

inline uint hash32(uint v) {
    v ^= v >> 16;
    v *= 0x7feb352du;
    v ^= v >> 15;
    v *= 0x846ca68bu;
    return v ^ (v >> 16);
}

inline float unit(uint seed, uint k) {
    return float(hash32(seed + 0x9e3779b9u * k) >> 8) * 0x1p-24f;
}

inline float log_sample(uint seed, uint k, float low, float high) {
    return exp2(mix(low, high, unit(seed, k)));
}

inline Point generate(uint index, uint salt) {
    const uint seed = hash32(index ^ (0x6ac690c5u * (salt + 1u)));
    const float offset = log_sample(seed, 1, -13.28771238f, 4.95419631f);
    const float m = 1.0f + offset;
    const float threshold = m * log(m) / (m - 1.0f);
    const float margin = log_sample(seed, 2, -13.28771238f, 2.0f);
    const float a = sqrt(threshold + margin);
    const float r = log_sample(seed, 3, -5.32192809f, 1.0f);
    const float vFraction = log_sample(seed, 4, -13.28771238f, -0.19264508f);
    const float r2 = r*r;
    const float c = vFraction*m / (r2*r2);
    const float eps = log_sample(seed, 5, -39.86313714f, -6.64385619f);
    const uint family = index & 7u;
    float offsetY = 0.0f;
    if (family >= 2u && family <= 4u) {
        offsetY = (2.0f * unit(seed, 6) - 1.0f) * r;
    } else if (family >= 5u && family <= 6u) {
        offsetY = (2.0f * unit(seed, 6) - 1.0f) * min(0.01f, 0.25f*r);
    } else if (family == 7u) {
        offsetY = 2.0f * (2.0f * unit(seed, 6) - 1.0f) * max(r, 0.5f);
    }
    return Point{m, a, c, r, eps, offsetY, margin};
}

// Integrate u on [origin+lo, origin+hi]. For a ridge panel origin is
// the candidate center and the signed displacement is computed without
// subtracting two nearly identical numbers.
inline float panel(Point p, constant float2* nodes, uint count,
                   float origin, float lo, float hi,
                   int side, uint power, float offset, float prefactor) {
    if (hi <= lo) return 0.0f;
    float total = 0.0f;
    const float width = hi - lo;
    const float direction = float(side);
    const float centeredDelta = p.delta - direction*origin;
    for (uint i = 0; i < count; ++i) {
        const float t = nodes[i].x;
        const float t2 = t*t;
        const float coordinate = power == 1 ? t : t2*t;
        const float jacobian = width * (power == 1 ? 1.0f : 3.0f*t2);
        const float s = lo + width*coordinate;
        const float u = origin + s;
        const float x = p.a + direction*u;
        const float dx = centeredDelta - direction*s;
        const float u2 = u*u;
        const float variance = p.epsilon + p.c*u2*u2;
        const float exponent = offset - 0.5f * (x*x + dx*dx/variance);
        total += nodes[i].y * jacobian * rsqrt(variance) * exp(exponent);
    }
    return prefactor * total;
}

inline float integrate(Point p, constant float2* nodes, uint count,
                       int side, float offset, float prefactor) {
    const float center = float(side) * p.delta;
    if (center <= 0.0f) {
        const float boundary = min(8.0f * sqrt(p.epsilon), p.r);
        return panel(p, nodes, count, 0, 0, boundary, side, 1, offset, prefactor)
             + panel(p, nodes, count, boundary, 0, p.r-boundary, side, 3, offset, prefactor);
    }
    const float center2 = center*center;
    const float scale = sqrt(p.epsilon + p.c*center2*center2);
    const float lo = max(-8.0f*scale, -center);
    const float hi = min(8.0f*scale, p.r-center);
    const float before = clamp(center+lo, 0.0f, p.r);
    const float after = clamp(center+hi, 0.0f, p.r);
    return panel(p, nodes, count, 0, 0, before, side, 3, offset, prefactor)
         + panel(p, nodes, count, center, lo, hi, side, 1, offset, prefactor)
         + panel(p, nodes, count, after, 0, p.r-after, side, 3, offset, prefactor);
}

kernel void quartic(constant float2* nodes [[buffer(0)]],
                    device Row* output [[buffer(1)]],
                    constant uint4& parameters [[buffer(2)]],
                    uint position [[thread_position_in_grid]],
                    uint lane [[thread_index_in_simdgroup]]) {
    const uint index = parameters.z + position;
    const Point p = generate(index, parameters.w);
    const float y = p.a + p.delta;
    const float variance = p.m + p.epsilon;
    const float offset = 0.5f*y*y/variance;
    const float prefactor = sqrt(variance * 0.15915494309189535f);
    const float ratio = integrate(p, nodes, parameters.y, -1, offset, prefactor)
                      + integrate(p, nodes, parameters.y, 1, offset, prefactor);
    const float threshold = p.a*p.a - p.margin;
    const bool yOutside = y*y > threshold;
    const bool atContact = p.delta == 0;
    const bool valid = isfinite(ratio);
    const uint group = position / 32u;
    for (uint stratum = 0; stratum < 8; ++stratum) {
        const bool enabled = stratum == 0
            || (stratum == 1 && yOutside)
            || (stratum == 2 && atContact)
            || (stratum == 3 && atContact && p.margin >= 0.5f)
            || (stratum == 4 && atContact && p.margin >= 1.0f)
            || (stratum == 5 && atContact && p.margin >= 0.5f && p.r <= 0.25f)
            || (stratum == 6 && atContact && p.margin >= 0.5f && p.r <= 0.125f)
            || (stratum == 7 && atContact && p.margin >= 0.5f && p.epsilon <= 1e-8f);
        const float value = (enabled && valid) ? ratio : -INFINITY;
        const float highest = simd_max(value);
        const uint winner = simd_min(value == highest ? lane : 32u);
        const uint tested = simd_sum(uint(enabled));
        const uint above = simd_sum(uint(enabled && valid && ratio > 1.0f));
        const uint nonfinite = simd_sum(uint(enabled && !valid));
        if (lane == winner) {
            output[group*8u+stratum] = Row{highest, p.m, p.a, p.c, p.r, p.epsilon,
                                          p.delta, p.margin, tested, above, nonfinite, index};
        }
    }
}
