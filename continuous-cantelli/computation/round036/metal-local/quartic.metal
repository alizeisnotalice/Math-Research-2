#include <metal_stdlib>
using namespace metal;

struct Point {
    float m, a, r, baseEdge, shape, power, epsilon, delta, margin, left, right;
    uint family;
};
struct Row {
    Point point;
    float ratio, baseline;
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
    const uint seed = hash32((index ^ 0xf6bacda7u) + 0x6ac690c5u*salt);
    const float m = 1.0f + log_sample(seed, 1, -13.28771238f, 4.95419631f);
    const float threshold = m*log(m)/(m-1.0f);
    const float margin = log_sample(seed, 2, -2.736965594f, 3.169925001f);
    const float a = sqrt(threshold+margin);
    const float r = log_sample(seed, 3, -5.643856190f, 0.0f);
    const uint family = index % 3u;
    const float left = mix(-0.875f,4.0f,unit(seed,4));
    const float right = family == 0u ? left : family == 1u ? 0.0f : mix(-0.875f,4.0f,unit(seed,5));
    const float q = log_sample(seed, 6, -5.0f, 6.0f);
    const float baseEdge = min(q*r*r, 0.8f*m/max(1.0f,1.0f+max(left,right)));
    const float shape = log_sample(seed, 7, -3.0f, 4.0f);
    const float power = index % 16u == 0u ? 2.0f : log_sample(seed, 8, -2.321928095f, 0.925999419f);
    const float epsilon = log_sample(seed, 9, -59.79470571f, -6.64385619f);
    const uint slot = (index / 3u) & 7u;
    float delta = 0.0f;
    if (slot >= 4u && slot <= 5u) delta = (2.0f*unit(seed,10)-1.0f)*0.04f*r;
    if (slot >= 6u) delta = (2.0f*unit(seed,10)-1.0f)*r;
    return Point{m,a,r,baseEdge,shape,power,epsilon,delta,margin,left,right,family};
}

inline float var(Point p, float u, float q) {
    const float t = clamp(u/p.r, 0.0f, 1.0f);
    const float t2 = t*t;
    float correction = 0;
    if (q != 0 && t > 0) {
        correction = q*exp(p.shape*(1.0f-powr(1.0f/t,p.power)));
    }
    return p.epsilon+p.baseEdge*t2*t2*(1.0f+correction);
}

inline float panel(Point p, constant float2* nodes, uint n, float lo, float hi,
                   float center, float q, float side, float source, bool linear) {
    if (hi <= lo) return 0;
    const float width = hi-lo;
    const bool forward = abs(lo) <= abs(hi);
    const float anchor = forward ? lo : hi;
    const float orientation = forward ? width : -width;
    float result = 0;
    for (uint k=0; k<n; ++k) {
        const float t = nodes[k].x;
        const float t2 = t*t;
        const float w = anchor + orientation*(linear ? t : t2*t);
        const float u = clamp(center+w,0.0f,p.r);
        const float x = p.a+side*u;
        const float inv = rsqrt(var(p,u,q));
        const float z = w*inv;
        result += nodes[k].y * (linear ? 1.0f : 3.0f*t2) * inv * exp(source-0.5f*(x*x+z*z));
    }
    return width*result;
}

inline float half_density(Point p, constant float2* nodes, uint n,
                          float side, float source, float q) {
    const float center = side*p.delta;
    const float lo = -center, hi = p.r-center;
    const float radius = 9.0f*sqrt(var(p,clamp(center,0.0f,p.r),q));
    const float ridgeLo = clamp(-radius,lo,hi), ridgeHi = clamp(radius,lo,hi);
    const float mid = p.r*p.r/sqrt(p.baseEdge);
    const float correction = q == 0.0f ? 0.0f : p.r*powr(1.0f+3.0f/p.shape,-1.0f/p.power);
    float cuts[8] = {lo, hi, ridgeLo, ridgeHi, clamp(0.0f,lo,hi),
                     clamp(0.3f*mid-center,lo,hi), clamp(mid-center,lo,hi),
                     clamp(correction-center,lo,hi)};
    for (uint j=1; j<8; ++j) {
        float value = cuts[j];
        uint k=j;
        while (k>0 && value<cuts[k-1]) {
            cuts[k] = cuts[k-1];
            --k;
        }
        cuts[k] = value;
    }
    float total = 0;
    for (uint j=1; j<8; ++j) {
        const bool linear = cuts[j]<=ridgeHi && cuts[j-1]>=ridgeLo;
        total += panel(p,nodes,n,cuts[j-1],cuts[j],center,q,side,source,linear);
    }
    return total;
}

kernel void quartic(constant float2* nodes [[buffer(0)]],
                    device Row* output [[buffer(1)]],
                    constant uint4& parameters [[buffer(2)]],
                    uint index [[thread_position_in_grid]],
                    uint lane [[thread_index_in_simdgroup]]) {
    const uint global = parameters.z+index;
    const Point p = generate(global,parameters.w);
    const float target = p.m+p.epsilon;
    const float y = p.a+p.delta;
    const float source = 0.5f*y*y/target;
    const float normal = sqrt(target * 0.15915494309189533577f);
    const float leftBase = half_density(p,nodes,parameters.y,-1.0f,source,0.0f);
    const float rightBase = half_density(p,nodes,parameters.y,1.0f,source,0.0f);
    const float baseline = normal*(leftBase+rightBase);
    const float ratio = normal*((p.left == 0 ? leftBase : half_density(p,nodes,parameters.y,-1.0f,source,p.left))
                 + (p.right == 0 ? rightBase : half_density(p,nodes,parameters.y,1.0f,source,p.right)));
    const bool valid = isfinite(ratio) && isfinite(baseline);
    const bool at = p.delta == 0;
    const bool upgrade = baseline <= 0.98f;
    const bool tight = at && p.r <= 0.125f && p.margin >= 0.5f && p.epsilon <= 1e-8f;
    const uint group = index/32u;
    for (uint stratum=0; stratum<8; ++stratum) {
        const bool enabled = stratum == 0
            || (stratum == 1 && at && upgrade)
            || (stratum == 2 && tight && upgrade && p.family == 0u)
            || (stratum == 3 && tight && upgrade && p.family != 0u)
            || (stratum == 4 && tight && upgrade && p.left >= 0 && p.right >= 0)
            || (stratum == 5 && tight && upgrade && p.left <= 0 && p.right <= 0)
            || (stratum == 6 && at && upgrade && p.m >= 1.5f && p.m <= 2.5f && p.a >= 1.4f)
            || (stratum == 7 && tight && upgrade && p.power == 2.0f);
        const float value = enabled && valid ? ratio : -INFINITY;
        const float best = simd_max(value);
        const uint winner = simd_min(best == value ? lane : 32u);
        const uint tested = simd_sum(uint(enabled));
        const uint above = simd_sum(uint(enabled && valid && ratio > 1.0f));
        const uint nonfinite = simd_sum(uint(enabled && !valid));
        if (lane == winner) {
            output[group*8u+stratum] = Row{p,best,baseline,tested,above,nonfinite,global,parameters.w};
        }
    }
}
