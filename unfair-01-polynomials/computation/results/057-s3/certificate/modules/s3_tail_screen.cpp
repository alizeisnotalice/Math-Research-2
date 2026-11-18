// Exact Cayley halfplane cuts and rigorously enclosed Schur descent for S={3}.
#include <algorithm>
#include <array>
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <gmpxx.h>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Int = std::int64_t;
using Poly = std::array<Int, 9>;
using Counts = std::array<std::uint64_t, 12>;

struct Group {
    Int L, T, N, D, b, c, dlo, dhi, elo, ehi, zlo, zhi, count;

    bool read(std::istream& stream) {
        if (!(stream >> L)) return false;
        if (!(stream >> T >> N >> D >> b >> c >> dlo >> dhi >> elo >> ehi >> zlo >> zhi >> count))
            throw std::runtime_error("truncated group");
        if (L < -4 || L > 4 || T < -104 || T > 104 || D < 5 || D > 468
            || N < -117 || N > 2702 || T*T-4*N != D || b < -120 || b > 120
            || c < -92 || c > 92 || dlo < -41 || dhi > 41 || elo < -9 || ehi > 9
            || dlo > dhi || elo > ehi || zlo > zhi || count <= 0 || count > 1577)
            throw std::runtime_error("invalid group values");
        return true;
    }
};

Int root_floor(Int n) {
    // The float supplies an initial estimate only; the inequalities below
    // determine the exact floor, including on a square.
    Int k = static_cast<Int>(std::sqrt(static_cast<double>(n)));
    while (k * k > n) --k;
    while ((k + 1) * (k + 1) <= n) ++k;
    return k;
}

Int choose(int n, int k) {
    if (k < 0 || k > n) return 0;
    Int binomial = 1;
    for (int j = 1; j <= k; ++j) binomial = binomial*(n-j+1)/j;
    return binomial;
}

constexpr Int power(Int x, int n) {
    Int out = 1;
    for (int i = 0; i < n; ++i) out *= x;
    return out;
}

// ceil(n/d), exact in signed arithmetic for positive denominator d.
Int ceil_div(Int n, Int d) {
    assert(d > 0);
    return n / d + (n % d > 0);
}

Int floor_div(Int n, Int d) {
    return -ceil_div(-n, d);
}

struct Halfspace {
    Int d, e, rhs;

    void clip(Int d_value, Int& lower, Int& upper) const {
        const Int hm = rhs - d*d_value;
        if (e > 0) lower = std::max(lower, ceil_div(hm, e));
        else if (e < 0) upper = std::min(upper, floor_div(-hm, -e));
        else if (hm > 0) lower = upper + 1;
    }
};

using Matrix = std::array<Poly, 9>;

Matrix cayley_matrix() {
    Matrix matrix{};
    for (int j = 0; j <= 8; ++j)
        for (int k = 0; k <= 8; ++k)
            for (int i = std::max(0, k - 8 + j); i <= std::min(j, k); ++i)
                matrix[j][k] += (Int(1) << j) * power(-1,k-i) * choose(j,i) * choose(8-j,k-i);
    return matrix;
}

std::vector<Halfspace> cuts(Group g, bool with_cayley) {
    const Poly doubled{2, 2*g.L, 2*g.L*g.L, g.T, 2*g.b, 2*g.c, 0, 0, 2};
    const std::array<std::pair<Int,Int>, 9> points{
        {{1,1}, {1,2}, {2,1}, {1,4}, {4,1}, {3,4}, {4,3}, {6,5}, {5,6}}
    };
    std::vector<Halfspace> result;
    result.reserve(with_cayley ? 27 : 9);
    for (auto [p, q] : points) {
        Poly weights{};
        for (int j = 0; j <= 8; ++j) weights[j] = power(p,j)*power(q,8-j);
        const Int base = std::inner_product(doubled.begin(),doubled.end(),weights.begin(),Int(0));
        const Int irr = weights[3];
        result.push_back({2*weights[6],2*weights[7],root_floor(irr*irr*g.D)+1-base});
    }
    if (!with_cayley) return result;
    const static Matrix matrix = cayley_matrix();
    for (int reverse = 0; reverse < 2; ++reverse) {
        const auto index = [&](int j) { return reverse ? 8-j : j; };
        for (int k = 0; k <= 8; ++k) {
            Int base = 0;
            for (int j = 0; j <= 8; ++j) base += doubled[j]*matrix[index(j)][k];
            const Int irr = matrix[index(3)][k];
            result.push_back({2*matrix[index(6)][k],2*matrix[index(7)][k],root_floor(irr*irr*g.D)+1-base});
        }
    }
    return result;
}

struct Interval {
    double lo, hi;
};

double down(double value) { return std::nextafter(value, -std::numeric_limits<double>::infinity()); }
double up(double value) { return std::nextafter(value, std::numeric_limits<double>::infinity()); }

Interval operator*(Interval x, Interval y) {
    const auto [lo, hi] = std::minmax({x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi});
    return {down(lo), up(hi)};
}

Interval operator-(Interval x, Interval y) { return {down(x.lo-y.hi),up(x.hi-y.lo)}; }

Interval scale_by_two(Interval x, int power_of_two) {
    return {down(std::scalbn(x.lo,power_of_two)), up(std::scalbn(x.hi,power_of_two))};
}

enum class Decision { pass, fail, fallback };

Decision interval_schur(std::array<Interval,9> asc) {
    for (int j = 8; j > 0; --j) {
        const Interval lead = asc[j];
        const Interval tail = asc[0];
        const double abs_low = std::max({0., tail.lo, -tail.hi});
        const double abs_high = std::max(-tail.lo,tail.hi);
        if (lead.hi <= abs_low) return Decision::fail;
        if (!(lead.lo > abs_high)) return Decision::fallback;
        std::array<Interval,9> next{};
        double magnitude = 0.;
        for (int i = 0; i < j; ++i) {
            next[i] = lead*asc[i+1] - tail*asc[j-1-i];
            magnitude = std::max({magnitude,std::abs(next[i].lo),std::abs(next[i].hi)});
        }
        if (!std::isfinite(magnitude) || magnitude == 0.) return Decision::fallback;
        const int exponent = std::ilogb(magnitude);
        for (int i = 0; i < j; ++i) asc[i] = scale_by_two(next[i], -exponent);
    }
    return Decision::pass;
}

// (u + v sqrt(D)) represents exact coefficients of consistently scaled
// real polynomials. conjugacy is selected only in the starting coefficient.
struct Quadratic {
    mpz_class u, v;
};

Quadratic times(const Quadratic& x, const Quadratic& y, long D) {
    return {x.u*y.u+D*x.v*y.v,x.u*y.v+x.v*y.u};
}

Quadratic minus(const Quadratic& x, const Quadratic& y) {
    return {x.u-y.u,x.v-y.v};
}

int sign(const Quadratic& x, long D) {
    const int su = mpz_sgn(x.u.get_mpz_t()), sv = mpz_sgn(x.v.get_mpz_t());
    if (sv == 0 || su == sv) return su;
    if (su == 0) return sv;
    const mpz_class square_difference = x.u*x.u-D*x.v*x.v;
    return su*mpz_sgn(square_difference.get_mpz_t());
}

bool exact_schur(const Poly& doubled, Group g, bool reverse, int embedding) {
    std::array<Quadratic,9> a{};
    for (int j = 0; j <= 8; ++j) {
        const int k = reverse ? 8-j : j;
        a[j] = {long(doubled[k]*(Int(1)<<j)), long((k==3 ? embedding : 0)*(Int(1)<<j))};
    }
    for (int degree = 8; degree > 0; --degree) {
        const Quadratic& lead = a[degree];
        const Quadratic& tail = a[0];
        if (sign(minus(lead,tail),g.D) <= 0 || sign({lead.u+tail.u,lead.v+tail.v},g.D) <= 0)
            return false;
        std::array<Quadratic,9> next;
        for (int j = 0; j < degree; ++j)
            next[j] = minus(times(lead,a[j+1],g.D),times(tail,a[degree-1-j],g.D));
        a = std::move(next);
    }
    return true;
}

Decision enclosure(const Poly& doubled, Group g, bool reverse, int embedding) {
    // Every endpoint is exactly representable: |scaled| is below 2^53.
    constexpr int bits = 23;
    const Int radical_floor = root_floor(g.D*(Int(1)<<(2*bits)));
    const double lo = std::ldexp(static_cast<double>(radical_floor),-bits);
    const double hi = std::ldexp(static_cast<double>(radical_floor+1),-bits);
    std::array<Interval,9> a{};
    for (int j = 0; j <= 8; ++j) {
        const int k = reverse ? 8-j : j;
        const double base = static_cast<double>(doubled[k]);
        if (k == 3)
            a[j] = embedding > 0 ? Interval{base+lo,base+hi} : Interval{base-hi,base-lo};
        else
            a[j] = {base,base};
        a[j].lo = std::ldexp(a[j].lo,j);
        a[j].hi = std::ldexp(a[j].hi,j);
    }
    return interval_schur(a);
}

bool square(Int value) {
    const Int k = root_floor(value);
    return k*k==value;
}

int first_failed(Group g, Int d, Int e, Counts& totals, bool exact_all) {
    const Poly doubled{2,2*g.L,2*g.L*g.L,g.T,2*g.b,2*g.c,2*d,2*e,2};
    for (int stage = 0; stage < 4; ++stage) {
        const bool reversed = stage >= 2;
        const int embedding = stage % 2 == 0 ? 1 : -1;
        const auto screened = enclosure(doubled,g,reversed,embedding);
        if (screened == Decision::fallback) ++totals[10];
        if (exact_all || screened == Decision::fallback) {
            ++totals[11];
            const bool pass = exact_schur(doubled,g,reversed,embedding);
            if (screened != Decision::fallback && (screened == Decision::pass) != pass)
                throw std::runtime_error("outward enclosure mismatched exact Schur");
            if (!pass) return stage;
        } else if (screened == Decision::fail) return stage;
    }
    return 4;
}

int main(int argc, char** argv) {
    static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53);
    if (std::fegetround() != FE_TONEAREST) throw std::runtime_error("unexpected floating-point rounding mode");
    if (argc < 2 || argc > 3 || (argc==3 && std::string(argv[2])!="--exact-all"))
        throw std::runtime_error("usage: s3_tail_screen GROUPS [--exact-all]");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot read group file");
    const bool exact_all = argc == 3;
    Counts total{};
    Group g;
    std::tuple<Int,Int,Int,Int,Int> previous{-5,0,0,0,0};
    while (g.read(input)) {
        if (square(g.D)) throw std::runtime_error("square discriminant");
        const auto key = std::tuple{g.L,g.T,g.D,g.b,g.c};
        if (!(previous < key)) throw std::runtime_error("unordered group or duplicate");
        previous = key;
        ++total[0];
        auto hh = cuts(g,true);
        Int enumerated = 0;
        for (Int d=g.dlo; d<=g.dhi; ++d) {
            Int lo = std::max(g.elo,g.zlo+g.L*d), hi = std::min(g.ehi,g.zhi+g.L*d);
            if (lo>hi) continue;
            enumerated += hi-lo+1;
            for (std::size_t k=0; k<9 && lo<=hi; ++k) hh[k].clip(d,lo,hi);
            if (lo>hi) continue;
            total[2] += hi-lo+1;
            for (std::size_t k=9; k<hh.size() && lo<=hi; ++k) hh[k].clip(d,lo,hi);
            if (lo>hi) continue;
            total[3] += hi-lo+1;
            for (Int e=lo; e<=hi; ++e) {
                const int stage = first_failed(g,d,e,total,exact_all);
                ++total[4+stage];
                if (stage == 4)
                    std::cout << '[' << g.L << ',' << g.L*g.L << ',' << g.T << ',' << g.N << ','
                              << g.D << ',' << g.b << ',' << g.c << ',' << d << ',' << e << "]\n";
            }
        }
        if (enumerated != g.count) throw std::runtime_error("Newton interval count mismatch");
        total[1] += enumerated;
    }
    if (!input.eof()) throw std::runtime_error("group input parse error");
    if (total[3]!=total[4]+total[5]+total[6]+total[7]+total[8])
        throw std::runtime_error("Schur counts mismatch");
    const std::array<const char*,12> labels{"groups","mothers","positive","cayley","outer_plus_reject",
        "outer_minus_reject","inner_plus_reject","inner_minus_reject","pass","reserved",
        "interval_fallback_calls","exact_calls"};
    std::cerr << "{";
    for (std::size_t k=0; k<labels.size(); ++k)
        std::cerr << (k==0 ? "" : ",") << '"' << labels[k] << "\":" << total[k];
    std::cerr << "}\n";
}
