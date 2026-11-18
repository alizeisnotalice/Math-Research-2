#include <metal_stdlib>
using namespace metal;

struct Point {
    float m, a, r, edge, shape, power, epsilon, delta, margin;
    uint family;
};
struct Row {
    Point point;
    float ratio;
    uint tested, above, nonfinite, index, salt;
};

inline uint hash32(uint v) {
    v ^= v >> 16;
    v *= 0x7feb352du;
    v ^= v >> 15;
    v *= 0x846ca68bu;
    return v ^ (v >> 16);
}

inline float unit(uint seed, uint k) {
    return float(hash32(seed + 0x9e3779b9u*k) >> 8) * 0x1p-24f;
}

inline float log_sample(uint seed, uint k, float lo, float hi) {
    return exp2(mix(lo, hi, unit(seed,k)));
}

inline Point generate(uint index, uint salt) {
    const uint seed = hash32(index ^ (0x6ac690c5u * (salt + 1u)));
    const float m = 1.0f + log_sample(seed, 1u, -13.28771238f, 4.95419631f);
    const float threshold = m*log(m)/(m-1.0f);
    const float margin = log_sample(seed, 2u, -13.28771238f, 2.0f);
    const float a = sqrt(threshold+margin);
    const float r = log_sample(seed, 3u, -6.0f, 1.0f);
    const float edge = min(0.875f*m, r*r*log_sample(seed, 4u, -9.0f, 8.0f));
    const float shape = log_sample(seed, 5u, -3.0f, 5.0f);
    const float power = log_sample(seed, 6u, -1.0f, 3.0f);
    const float epsilon = log_sample(seed, 7u, -66.4385619f, -6.6438562f);
    const uint family = index % 3u;
    const uint slot = (index / 3u) & 7u;
    float delta = 0.0f;
    if (slot >= 3u && slot <= 5u) {
        delta = (2.0f*unit(seed,8u)-1.0f)*r;
    } else if (slot >= 6u) {
        delta = (2.0f*unit(seed,8u)-1.0f)*0.04f*r;
    }
    return Point{m,a,r,edge,shape,power,epsilon,delta,margin,family};
}

inline float variance(Point p, float u, float shape) {
    if (shape == 0.0f || u <= 0.0f) return p.epsilon;
    const float ratio = max(p.r/u, 1.0f);
    const float power = powr(ratio, p.power);
    return p.epsilon + p.edge * exp(-shape*(power-1.0f));
}

inline float variance_level(Point p, float drop, float shape) {
    return p.r * powr(1.0f + max(0.0f, drop)/shape, -1.0f/p.power);
}

inline float panel(Point p, constant float2* nodes, uint n,
                 float lo, float hi, float center, float shape,
                 float side, float sourceOffset, bool linear) {
    if (hi <= lo) return 0.0f;
    const float width = hi-lo;
    const bool startLow = abs(lo) <= abs(hi);
    const float anchor = startLow ? lo : hi;
    const float orientation = startLow ? width : -width;
    float sum = 0.0f;
    for (uint k=0; k<n; ++k) {
        const float t = nodes[k].x;
        const float t2 = t*t;
        const float w = anchor + orientation * (linear ? t : t2*t);
        const float u = clamp(center+w, 0.0f, p.r);
        const float x = p.a + side*u;
        const float inv = rsqrt(variance(p, u, shape));
        const float z = w*inv;
        const float ex = sourceOffset - 0.5f*(x*x+z*z);
        const float jac = width * (linear ? 1.0f : 3.0f*t2);
        sum += nodes[k].y*jac*inv*exp(ex);
    }
    return sum;
}

inline float half_integral(Point p, constant float2* nodes, uint n,
                  float side, float sourceOffset) {
    const float shape = (p.family == 1u && side > 0.0f) ? 0.0f
                       : (p.family == 2u && side > 0.0f) ? 2.0f*p.shape : p.shape;
    const float center = side*p.delta;
    const float lo = -center, hi = p.r-center;
    const float ridge = 9.0f*sqrt(variance(p, clamp(center, 0.0f, p.r), shape));
    const float ridgeLo = clamp(-ridge, lo, hi), ridgeHi = clamp(ridge, lo, hi);
    if (shape == 0.0f) {
        return panel(p,nodes,n,ridgeLo,ridgeHi,center,shape,side,sourceOffset,true);
    }
    const float drop = log(p.edge/(p.r*p.r));
    float x[8] = {lo, hi, ridgeLo, ridgeHi, clamp(0.0f,lo,hi),
                  clamp(variance_level(p,drop+12.0f,shape)-center,lo,hi),
                  clamp(variance_level(p,drop+4.0f,shape)-center,lo,hi),
                  clamp(variance_level(p,drop-4.0f,shape)-center,lo,hi)};
    for (uint j=1; j<8; ++j) {
        float value = x[j];
        uint k=j;
        while (k > 0 && value < x[k-1]) {
            x[k] = x[k-1];
            --k;
        }
        x[k] = value;
    }
    float total = 0.0f;
    for (uint j=1; j<8; ++j) {
        const bool linear = (x[j] <= ridgeHi && x[j-1] >= ridgeLo);
        total += panel(p,nodes,n,x[j-1],x[j],center,shape,side,sourceOffset,linear);
    }
    return total;
}

kernel void flat_core(constant float2* nodes [[buffer(0)]],
                      device Row* output [[buffer(1)]],
                      constant uint4& arguments [[buffer(2)]],
                      uint index [[thread_position_in_grid]],
                      uint lane [[thread_index_in_simdgroup]]) {
    const uint global = arguments.z + index;
    const Point p = generate(global, arguments.w);
    const float target = p.m+p.epsilon;
    const float y = p.a+p.delta;
    const float sourceOffset = 0.5f*y*y/target;
    const float ratio = sqrt(target * 0.15915494309189533577f) *
        (half_integral(p,nodes,arguments.y,-1.0f,sourceOffset)+half_integral(p,nodes,arguments.y,1.0f,sourceOffset));
    const bool valid = isfinite(ratio);
    const bool at = p.delta == 0.0f;
    const bool tight = at && p.margin >= 0.25f && p.r <= 0.125f && p.epsilon <= 1e-8f;
    const uint group = index/32u;
    for (uint stratum=0; stratum<8; ++stratum) {
        const bool enabled = stratum == 0
            || (stratum == 1 && y*y >= p.a*p.a-p.margin)
            || (stratum == 2 && at && p.family == 0u && p.margin >= 0.5f)
            || (stratum == 3 && at && p.family == 1u && p.margin >= 0.25f)
            || (stratum == 4 && tight && p.family == 0u)
            || (stratum == 5 && tight && p.family == 1u)
            || (stratum == 6 && tight && p.family == 2u)
            || (stratum == 7 && at && p.m >= 1.5f && p.m <= 2.5f && p.a >= 1.4f && p.r <= 0.125f);
        const float value = enabled && valid ? ratio : -INFINITY;
        const float highest = simd_max(value);
        const uint winner = simd_min(value == highest ? lane : 32u);
        const uint tested = simd_sum(uint(enabled));
        const uint above = simd_sum(uint(enabled && valid && ratio>1.0f));
        const uint nonfinite = simd_sum(uint(enabled && !valid));
        if (lane == winner) {
            output[group*8u+stratum] = Row{p,highest,tested,above,nonfinite,global,arguments.w};
        }
    }
}
