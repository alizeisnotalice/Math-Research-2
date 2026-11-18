// Strict positivity of both real embeddings at rational points.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;
struct Group {
    Int L, M, T, N, D, b, c, dlo, dhi, elo, ehi, zlo, zhi, count;
};
struct Point {
    Int p, q;
    std::array<Int, 9> weights;
    std::array<Int, 469> threshold;
};
struct Plane { Int d_weight, e_weight, rhs; };

Int power(Int x, int n) {
    Int value = 1;
    while (n--) value *= x;
    return value;
}

Int integer_sqrt(Int n) {
    Int low = 0, high = 1;
    while (Wide(high)*high <= n) high *= 2;
    while (high-low > 1) {
        Int mid = (low+high)/2;
        if (Wide(mid)*mid <= n) low = mid;
        else high = mid;
    }
    return low;
}

Point point(Int p, Int q) {
    Point r{p, q, {}, {}};
    for (int i = 0; i <= 8; ++i) r.weights[i] = power(p,i)*power(q,8-i);
    for (int D = 5; D < static_cast<int>(r.threshold.size()); ++D) {
        const Int w = r.weights[3];
        r.threshold[D] = integer_sqrt(w*w*D)+1;
    }
    return r;
}

Int amount(const Group& g, const Point& t) {
    const auto& w = t.weights;
    return 2*(w[0]+w[1]*g.L+w[2]*g.M+w[4]*g.b+w[5]*g.c+w[8])+g.T*w[3];
}

Int div_ceil(Int a, Int b) {
    const Int quotient = a/b;
    return quotient + (a%b > 0);
}

std::vector<Plane> planes(const Group& g, const std::vector<Point>& points) {
    std::vector<Plane> result;
    for (const Point& t : points) {
        result.push_back({2*t.weights[6], 2*t.weights[7], t.threshold[g.D]-amount(g,t)});
    }
    return result;
}

std::pair<Int, Int> initial_e(const Group& g, Int d) {
    return {std::max(g.elo,g.zlo+g.L*d), std::min(g.ehi,g.zhi+g.L*d)};
}

Int clip(Int lower, Int d, const std::vector<Plane>& constraints) {
    for (const Plane& z : constraints) lower = std::max(lower,div_ceil(z.rhs-z.d_weight*d,z.e_weight));
    return lower;
}

bool direct(const Group& g, Int d, Int e, const Point& t) {
    // Separate Horner evaluation of twice the rational part and the radical weight.
    const std::array<Int,9> coeff{2,2*g.L,2*g.M,g.T,2*g.b,2*g.c,2*d,2*e,2};
    Int left = coeff[8], q_power = 1;
    for (int k = 7; k >= 0; --k) {
        q_power *= t.q;
        left = left * t.p + coeff[k] * q_power;
    }
    const Int radical = power(t.p,3)*power(t.q,5);
    return left > 0 && Wide(left)*left > Wide(g.D)*radical*radical;
}

bool read_group(Group& g) {
    if (!(std::cin >> g.L)) return false;
    if (!(std::cin >> g.M >> g.T >> g.N >> g.D >> g.b >> g.c >> g.dlo >> g.dhi >> g.elo >> g.ehi >> g.zlo >> g.zhi >> g.count))
        throw std::runtime_error("short or malformed tail group");
    if (g.D != g.T*g.T-4*g.N || g.D < 5 || g.D > 468 || g.M != g.L*g.L || g.L < -4 || g.L > 4 ||
        g.dlo < -41 || g.dlo > g.dhi || g.dhi > 41 || g.elo < -9 || g.elo > g.ehi || g.ehi > 9 ||
        g.zlo > g.zhi || g.count <= 0 || g.count > 1577)
        throw std::runtime_error("invalid group envelope");
    return true;
}

int main(int argc, char** argv) {
    if (argc > 2) throw std::invalid_argument("usage: s3_positive_axis [crosscheck_stride] < groups.txt");
    const Int stride = argc == 2 ? std::stoll(argv[1]) : 257;
    if (stride < 1) throw std::invalid_argument("positive crosscheck stride required");
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    const std::array<std::pair<Int,Int>,9> axis{{{1,1},{1,2},{2,1},{1,4},{4,1},{3,4},{4,3},{6,5},{5,6}}};
    std::vector<Point> points;
    for (auto [p,q] : axis) points.push_back(point(p,q));
    Int groups=0, retained_groups=0, mothers=0, pass=0, cross=0, exhausted_rows=0;
    Group g{};
    while (read_group(g)) {
        const auto constraints = planes(g,points);
        Int counted=0, accepted=0;
        for (Int d=g.dlo; d<=g.dhi; ++d) {
            const auto [lo,hi] = initial_e(g,d);
            if (lo>hi) continue;
            const Int bound = clip(lo,d,constraints);
            if (groups % stride == 0) {
                for (Int e=lo; e<=hi; ++e) {
                    const bool feasible = std::all_of(points.begin(),points.end(),[&](const Point& t){return direct(g,d,e,t);});
                    if (feasible != (e>=bound)) throw std::runtime_error("direct positivity disagreement");
                    ++cross;
                }
            }
            counted += hi-lo+1;
            if (bound<=hi) accepted+=hi-bound+1;
            else ++exhausted_rows;
        }
        if (counted != g.count) throw std::runtime_error("tail group count disagreement");
        ++groups;
        if (accepted) ++retained_groups;
        mothers += counted;
        pass += accepted;
    }
    if (!std::cin.eof() || groups==0) throw std::runtime_error("missing or malformed input");
    std::cout << "{\"groups\":" << groups << ",\"mothers\":" << mothers
              << ",\"axis_pass\":" << pass << ",\"axis_reject\":" << mothers-pass
              << ",\"retained_groups\":" << retained_groups << ",\"exhausted_d_rows\":" << exhausted_rows
              << ",\"direct_integer_crosschecks\":" << cross << "}\n";
}
