#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using I64 = std::int64_t;
using I128 = __int128_t;
using Box = std::array<I64, 11>;
constexpr I64 grid = I64(1) << 44;

struct RadialLine {
    I64 radius_num, radius_den, lower, center_upper, old_clip;

    I64 at(I64 distance) const {
        const I128 excess = I128(distance) * radius_den - I128(radius_num) * grid;
        if (excess < 0)
            return 0;
        return lower + I64(excess * (lower - center_upper) / (I128(radius_num) * grid));
    }
};

I64 sqrt_floor(I128 squared) {
    I64 distance = I64(std::sqrt(static_cast<long double>(squared)));
    while (I128(distance) * distance > squared)
        --distance;
    while (I128(distance + 1) * (distance + 1) <= squared)
        ++distance;
    return distance;
}

I64 sqrt_ceil(I128 squared) {
    const I64 f = sqrt_floor(squared);
    return f + (I128(f) * f < squared);
}

I128 sq(I64 x) { return I128(x) * x; }

struct Traversal {
    const std::uint8_t* codes;
    const std::vector<int>& end;
    const std::vector<int>& prefix;
    const RadialLine* lines;
    I64* scores;
    I64* exterior;

    bool visit(int node, Box& lo, Box& hi, int depth) const {
        if (depth > 256)
            return false;
        if (prefix[node] == prefix[end[node]])
            return true;
        const int axis = codes[node];
        if (axis >= 253) {
            const auto distance_to_zero = [&](int j) { return std::max({I64(0), lo[j], -hi[j]}); };
            I64 score = 0;
            I64 cap = std::numeric_limits<I64>::max();
            bool empty_clip = false;
            for (int b = 0; b < 4; ++b) {
                const int j = 3 * b;
                const I64 x = distance_to_zero(j), y = distance_to_zero(j + 1);
                const I128 closest = sq(x) + sq(y);
                score = std::max(score, lines[b].at(sqrt_floor(closest)));

                const I128 radius = sq(lines[b].old_clip);
                if (closest > radius) {
                    empty_clip = true;
                    continue;
                }
                const I64 x_cap = sqrt_ceil(radius - sq(y));
                const I64 y_cap = sqrt_ceil(radius - sq(x));
                if (lo[j] < -x_cap || hi[j] > x_cap)
                    cap = std::min(cap, lines[b].at(sqrt_floor(sq(x_cap) + sq(y))));
                if (lo[j + 1] < -y_cap || hi[j + 1] > y_cap)
                    cap = std::min(cap, lines[b].at(sqrt_floor(sq(y_cap) + sq(x))));
            }
            scores[prefix[node]] = score;
            exterior[prefix[node]] = empty_clip ? score : cap;
            return true;
        }

        const I64 lower = lo[axis], upper = hi[axis];
        if ((lower + upper) % 2 || upper - lower < 2)
            return false;
        const I64 middle = (lower + upper) / 2;
        hi[axis] = middle;
        if (!visit(node + 1, lo, hi, depth + 1))
            return false;
        lo[axis] = middle;
        hi[axis] = upper;
        if (!visit(end[node + 1], lo, hi, depth + 1))
            return false;
        lo[axis] = lower;
        return true;
    }
};

extern "C" int radial_scores(const std::uint8_t* codes, int nodes, const I64* low, const I64* high,
                             const std::uint8_t* active, const I64* params, I64* scores, I64* exterior) {
    if (nodes < 1 || nodes >= (1 << 28))
        return 1;
    std::vector<int> end(nodes), prefix(nodes + 1);
    for (int i = nodes - 1; i >= 0; --i) {
        if (codes[i] >= 253)
            end[i] = i + 1;
        else {
            if (codes[i] >= 11 || i + 1 >= nodes || end[i + 1] >= nodes)
                return 2;
            end[i] = end[end[i + 1]];
        }
    }
    if (end[0] != nodes)
        return 3;
    for (int i = 0; i < nodes; ++i) {
        if (active[i] > 1 || (active[i] && codes[i] < 253))
            return 4;
        prefix[i + 1] = prefix[i] + active[i];
    }
    Box lo, hi;
    std::copy(low, low + 11, lo.begin());
    std::copy(high, high + 11, hi.begin());
    RadialLine lines[4];
    for (int b = 0; b < 4; ++b)
        lines[b] = {params[5*b], params[5*b+1], params[5*b+2], params[5*b+3], params[5*b+4]};
    return Traversal{codes, end, prefix, lines, scores, exterior}.visit(0, lo, hi, 0) ? 0 : 5;
}
