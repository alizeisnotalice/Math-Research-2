// Exact Cayley halfplane cuts and rigorously enclosed Schur descent for S={3}.
#include <algorithm>
#include <array>
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <bit>
#include <stdexcept>
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

// All endpoints are exact signed integers in a common dyadic scale.
using Wide = __int128_t;
using UWide = __uint128_t;
struct IntBox { Int lo, hi; };
struct WideBox { Wide lo, hi; };

WideBox product(IntBox x, IntBox y) {
    auto [lo,hi] = std::minmax({Wide{x.lo}*y.lo,Wide{x.lo}*y.hi,Wide{x.hi}*y.lo,Wide{x.hi}*y.hi});
    return {lo,hi};
}

int width(UWide n) {
    const auto high = static_cast<std::uint64_t>(n >> 64);
    return high ? 64+std::bit_width(high) : std::bit_width(static_cast<std::uint64_t>(n));
}

int schur(std::array<IntBox,9> a) {
    for (int n=8;n>=1;--n) {
        const IntBox lead=a[n],tail=a[0];
        const Int max_abs=std::max(tail.hi,-tail.lo);
        const Int min_abs=std::max({Int{0},tail.lo,-tail.hi});
        if (lead.hi<=min_abs) return 0;
        if (lead.lo<=max_abs) return -1;
        if (n==1) return 1;
        std::array<WideBox,8> next{};
        Wide maximum=0;
        for (int j=0;j<n;++j) {
            const auto p=product(lead,a[j+1]), q=product(tail,a[n-1-j]);
            next[j]={p.lo-q.hi,p.hi-q.lo};
            maximum=std::max({maximum,next[j].hi,-next[j].lo});
        }
        if (!maximum) throw std::runtime_error("zero polynomial after strict step");
        const int shift=std::max(0,width(static_cast<UWide>(maximum))-51);
        for (int j=0;j<n;++j) {
            // C++20 signed shifts divide by 2^shift and round toward -infinity.
            a[j]={static_cast<Int>(next[j].lo>>shift),static_cast<Int>(-((-next[j].hi)>>shift))};
        }
    }
    throw std::runtime_error("unreachable degree");
}

constexpr int precision=40;
std::array<Int,469> root_table() {
    std::array<Int,469> table{};
    for (int d=5;d<=468;++d) {
        const Wide target=Wide{d}<<(2*precision);
        Int lo=0,hi=Int{32}<<precision;
        while (hi-lo>1) {
            const Int middle=lo+(hi-lo)/2;
            if (Wide{middle}*middle<=target) lo=middle;
            else hi=middle;
        }
        table[d]=lo;
    }
    return table;
}

int first_failed(const Group& g, Int d, Int e, Counts& ct) {
    const Poly base{2,2*g.L,2*g.L*g.L,g.T,2*g.b,2*g.c,2*d,2*e,2};
    static const auto radicals=root_table();
    const Int q=radicals[g.D];
    for (int stage=0;stage<4;++stage) {
        std::array<IntBox,9> x{};
        for (int j=0;j<=8;++j) {
            const int k=stage>=2?8-j:j;
            const Int rational=base[k]*(Int{1}<<precision);
            Int lo=rational,hi=rational;
            if (k==3) {
                if (stage%2) {lo-=q+1;hi-=q;}
                else {lo+=q;hi+=q+1;}
            }
            x[j]={lo*(Int{1}<<j),hi*(Int{1}<<j)};
        }
        const int r=schur(x);
        if (r==0) return stage;
        if (r==-1) {++ct[10];return 5;}
    }
    return 4;
}

bool square(Int value) {
    const Int k=root_floor(value);
    return k*k==value;
}

int main(int argc, char** argv) {
    if (argc != 2) throw std::runtime_error("usage: s3_integer_enclosure GROUPS");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot read group file");
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
                const int stage = first_failed(g,d,e,total);
                ++total[4+stage];
                if (stage >= 4)
                    std::cout << '[' << g.L << ',' << g.L*g.L << ',' << g.T << ',' << g.N << ','
                              << g.D << ',' << g.b << ',' << g.c << ',' << d << ',' << e << "]\n";
            }
        }
        if (enumerated != g.count) throw std::runtime_error("Newton interval count mismatch");
        total[1] += enumerated;
    }
    if (!input.eof()) throw std::runtime_error("group input parse error");
    if (total[3]!=total[4]+total[5]+total[6]+total[7]+total[8]+total[9])
        throw std::runtime_error("Schur counts mismatch");
    const std::array<const char*,12> labels{"groups","mothers","positive","cayley","outer_plus_reject",
        "outer_minus_reject","inner_plus_reject","inner_minus_reject","pass","unresolved",
        "interval_fallback_calls","exact_calls"};
    std::cerr << "{";
    for (std::size_t k=0; k<labels.size(); ++k)
        std::cerr << (k==0 ? "" : ",") << '"' << labels[k] << "\":" << total[k];
    std::cerr << "}\n";
}
