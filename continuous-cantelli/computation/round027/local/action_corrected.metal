#include <metal_stdlib>
using namespace metal;

struct Input {
    float logEta, r, hR, bR2, a, logLambda, logK;
    uint n;
};
struct Answer {
    float metric, t, realGap, firstGap, analyticGap, sourceGap, remainder, logAmplitude;
};
struct Row {
    Input input;
    Answer answer;
    uint tested, realAbove, correctedAbove, firstDisagrees, nonfinite, index, salt, padding;
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

inline float integer_power(float x, uint n) {
    float result = 1;
    while (n != 0) {
        if ((n & 1u) != 0) result *= x;
        x *= x;
        n >>= 1;
    }
    return result;
}

inline float power_delta(float t, uint n) {
    float value = 0;
    const float x = 1.0f+t;
    for (uint i = 0; i < n; ++i) value = x*value+t;
    return value;
}

inline float stable_log1p(float x) {
    if (abs(x) >= 0.04f) return log(1.0f+x);
    return x * (1.0f+x*(-0.5f+x*(1.0f/3.0f+x*(-0.25f+x*(0.2f-x/6.0f)))));
}

inline float stable_expm1(float x) {
    if (abs(x) >= 0.125f) return exp(x)-1.0f;
    return x * (1.0f+x*(0.5f+x*(1.0f/6.0f+x*(1.0f/24.0f+x*(1.0f/120.0f+x/720.0f)))));
}

inline Input generate(uint index, uint salt) {
    const uint seed = hash32(index ^ (0x6ac690c5u * (salt + 1u)));
    const uint n = 4u + 2u * (index % 3u);
    const uint family = (index / 3u) & 7u;
    const float logEta = mix(-16.0f, -0.5f, unit(seed, 1));
    const float logR = mix(-5.0f, -1.0f, unit(seed, 2));
    const float r = exp2(logR);
    const float hR = family == 0u || family == 1u ? 0.0f : 0.6f * (unit(seed, 3) - 0.5f);
    const float bR2 = 0.3f * (unit(seed, 4) - 0.5f);
    const float a = 12.0f * (unit(seed, 5) - 0.5f);
    // Enforce lambda>=64 and c*r^n<=1/2 using c=2/(n*lambda*s^(n-2)).
    const float logLambda = max(6.0f, 2.0f - log2(float(n)) + 2.0f * logR
                                  - float(n - 2u) * logEta) + 8.0f * unit(seed, 6);
    const float deltaK = exp2(mix(-23.0f, -1.5f, unit(seed, 7)));
    const float logK = family <= 5u ? 0.0f : stable_log1p(family == 6u ? deltaK : -deltaK);
    return Input{logEta, r, hR, bR2, a, logLambda, logK, n};
}

inline float scaled_derivative(float t, float U, float V, float H, float J, uint n) {
    const float x = 1.0f+t;
    const float nf = float(n);
    const float base = power_delta(t, 2u*(n-1u));
    return base + (1.0f+base)*x*((nf+1.0f)*U/nf + (nf+2.0f)*V*x/nf)
                + (H+J*x)*integer_power(x,n-1u);
}

inline float root(float U, float V, float H, float J, float xmax, uint n) {
    float lo = -0.75f;
    float hi = min(1.0f, xmax - 1.0f);
    const float center = scaled_derivative(0.0f, U, V, H, J, n);
    if (center == 0.0f) return 0.0f;
    const float span = min(1.0f, 2.0f*abs(center));
    if (center > 0) {
        hi = 0;
        const float candidate = max(lo, -span);
        if (scaled_derivative(candidate, U, V, H, J, n) <= 0.0f) lo = candidate;
    } else {
        lo = 0;
        const float candidate = min(hi, span);
        if (scaled_derivative(candidate, U, V, H, J, n) >= 0.0f) hi = candidate;
    }
    for (uint i = 0; i < 30u; ++i) {
        const float t = (lo+hi)*0.5f;
        const float value = scaled_derivative(t, U, V, H, J, n);
        if (value == 0.0f) return t;
        if (value > 0.0f) hi = t; else lo = t;
    }
    return (lo+hi)*0.5f;
}

// ((1+t)^n - 1 - n*t)/t^2, continued to t=0.
inline float second_polynomial(float x, uint n) {
    float p = 0.0f;
    for (uint j = 1; j < n; ++j) p = float(j) + x*p;
    return p;
}

inline float remainder(float t, float U, float V, uint n) {
    const uint m = n - 2u;
    const float x = 1.0f + t;
    const float mf = float(m);
    const float nf = float(n);
    const float baseline = t*t * (second_polynomial(x, n)/nf
        + (mf*mf + (mf*t-1.0f) * second_polynomial(x, m))/(mf*integer_power(x, m)));
    const float linear = U * power_delta(t, n+1u)/nf;
    return baseline + linear + V*integer_power(x, n+2u)/nf;
}

inline Answer evaluate(Input p) {
    const uint n = p.n;
    const float nf = float(n);
    const float theta = 1.57079632679489661923f/(nf-1.0f);
    const float A = (1.0f - 1.0f/nf)*sin(theta);
    const float q = nf/(2.0f*(nf-1.0f));
    const float baseGap = A * stable_expm1(q*p.logK);
    const float zPower = (nf-2.0f)*sin(theta)*0.5f*exp(q*p.logK);
    const float z0 = powr(zPower, 1.0f/nf);
    const float eta = exp2(p.logEta);
    const float etaZ = eta*z0;
    const float U = -abs(p.hR)*etaZ;
    const float V = p.bR2*etaZ*etaZ;
    const float xmax = 1.0f/etaZ;
    const float t = root(U, V, 0.0f, 0.0f, xmax, n);
    const float rem = remainder(t, U, V, n);
    const float realGap = baseGap + zPower*(U/nf+rem);
    const float linearOfAnalytic = -abs(p.hR)*eta*sin(2.0f*theta)/nf;
    const float quadraticOfAnalytic = -sin(3.0f*theta)*eta*eta*(p.bR2/nf
        - p.hR*p.hR*(nf+1.0f)*(nf+1.0f)/(2.0f*nf*nf*(nf-1.0f)));
    const float linearOfReal = zPower*U/nf;
    const float firstGap = baseGap + linearOfReal - linearOfAnalytic;
    const float analyticGap = firstGap + zPower*rem - quadraticOfAnalytic;

    const float u0 = etaZ*p.r;
    const float invLambda = exp2(-p.logLambda);
    const float H = p.a*u0*invLambda/zPower;
    const float J = u0*u0*invLambda/zPower;
    const float UH = p.hR * etaZ;
    const float tPlus = root(UH, V, H, J, xmax, n);
    const float tMinus = root(-UH, V, -H, J, xmax, n);
    const float pX = 1.0f+tPlus;
    const float mX = 1.0f+tMinus;
    const float plus = UH/nf+remainder(tPlus,UH,V,n)+H*pX+0.5f*J*pX*pX;
    const float minus = -UH/nf+remainder(tMinus,-UH,V,n)-H*mX+0.5f*J*mX*mX;
    const float sourceGap = baseGap + zPower*min(plus,minus);
    const float secondary = V/nf - ((nf+1.0f)*(nf+1.0f)/(4.0f*nf*nf*(nf-1.0f)))*U*U;
    const float metric = abs(rem-secondary) / max(1e-25f, abs(V)+U*U);
    return Answer{metric, t, realGap, firstGap, analyticGap, sourceGap, rem,
                  log(A+realGap)+p.logLambda*0.6931471805599453f};
}

kernel void action(device Row* output [[buffer(0)]],
                   constant uint2& parameters [[buffer(1)]],
                   uint position [[thread_position_in_grid]],
                   uint lane [[thread_index_in_simdgroup]]) {
    const uint index = parameters.x + position;
    const Input p = generate(index, parameters.y);
    const Answer a = evaluate(p);
    const bool critical = p.logK == 0.0f;
    const bool small = p.logEta <= -9.0f;
    const bool valid = isfinite(a.metric) && isfinite(a.realGap) && isfinite(a.sourceGap);
    const uint group = position / 32u;
    for (uint stratum = 0; stratum < 12; ++stratum) {
        const bool orderMatch = p.n == 4u + 2u * (stratum / 4u);
        const uint family = stratum % 4u;
        const bool enabled = orderMatch && ((family == 0u && critical && small && abs(p.hR) >= 0.125f)
                           || (family == 1u && critical && small && p.hR == 0.0f && abs(p.bR2) >= 0.0625f)
                           || (family == 2u && critical && !small)
                           || (family == 3u && !critical));
        const float metric = (enabled && valid) ? a.metric : -INFINITY;
        const float highest = simd_max(metric);
        const uint winner = simd_min(metric == highest ? lane : 32u);
        const uint tested = simd_sum(uint(enabled));
        const uint realAbove = simd_sum(uint(enabled && valid && a.realGap > 0.0f));
        const uint correctedAbove = simd_sum(uint(enabled && valid && a.analyticGap > 0.0f));
        const uint firstDisagrees = simd_sum(uint(enabled && valid && ((a.firstGap > 0.0f) != (a.analyticGap > 0.0f))));
        const uint nonfinite = simd_sum(uint(enabled && !valid));
        if (lane == winner) {
            Answer out = a;
            out.metric = highest;
            output[group*12u+stratum] = Row{p, out, tested, realAbove, correctedAbove,
                                            firstDisagrees, nonfinite, index, parameters.y, 0};
        }
    }
}
