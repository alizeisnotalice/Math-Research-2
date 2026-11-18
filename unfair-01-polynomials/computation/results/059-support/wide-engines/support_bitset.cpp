// Exact support-pair classification using shifted uint64_t bit sets.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Unsigned 128-bit exponent storage; all count arithmetic is still decimal JSON.
int trailing_zero(__uint128_t value) {
    const auto low = static_cast<std::uint64_t>(value);
    return low ? std::countr_zero(low) : 64+std::countr_zero(static_cast<std::uint64_t>(value>>64));
}

std::string decimal(__uint128_t value) {
    std::string digits;
    do {
        digits.push_back(static_cast<char>('0'+value%10));
        value/=10;
    } while (value);
    std::reverse(digits.begin(),digits.end());
    return digits;
}

using Mask = __uint128_t;
struct Support {
    Mask points;
    Mask reflected;
    int terms;
    int gcd;
};

void enumerate(std::vector<Support>& output, int degree, int left, int next,
               Mask points, Mask reflected, int gcd, int terms) {
    if (left == 0) {
        output.push_back({points, reflected, terms, gcd});
        return;
    }
    for (int i = next; i + left <= degree; ++i) {
        enumerate(output, degree, left - 1, i + 1, points | (Mask{1} << i),
                  reflected | (Mask{1} << (degree - i)), std::gcd(gcd, i), terms);
    }
}

std::vector<Support> supports(int degree, int cap) {
    std::vector<Support> output;
    const Mask endpoints = Mask{1} | (Mask{1} << degree);
    for (int k = 2; k <= std::min(degree + 1, cap); ++k) {
        enumerate(output, degree, k - 2, 1, endpoints, endpoints, degree, k);
    }
    return output;
}

bool before(int ac, Mask a, int bc, Mask b) {
    if (ac != bc) return ac < bc;
    const Mask diff = a ^ b;
    return diff && (a & (diff & -diff));
}

bool canonical(const Support& a, const Support& b, bool equal_degree) {
    if (equal_degree && before(b.terms, b.points, a.terms, a.points)) return false;
    if (std::gcd(a.gcd, b.gcd) != 1) return false;
    Mask ar = a.reflected, br = b.reflected;
    int ac = a.terms, bc = b.terms;
    if (equal_degree && before(bc, br, ac, ar)) {
        std::swap(ac, bc);
        std::swap(ar, br);
    }
    if (a.points == ar) return !before(bc, br, b.terms, b.points);
    return !before(ac, ar, a.terms, a.points);
}

// 0 excluded, 1 unique representations, 2 unresolved by this rule.
int classify(Mask a, Mask b, int m, int n) {
    Mask seen = 0, repeated = 0;
    for (Mask t = a; t; t &= t - 1) {
        const Mask shifted = b << trailing_zero(t);
        repeated |= seen & shifted;
        seen |= shifted;
    }
    if (!repeated) return 1;
    const Mask single = seen & ~repeated;
    Mask ua = Mask{1} | (Mask{1} << m);
    Mask ub = Mask{1} | (Mask{1} << n);
    for (Mask t = a; t; t &= t - 1) {
        const int i = trailing_zero(t);
        const Mask unique_partners = (single >> i) & b;
        if (unique_partners) ua |= Mask{1} << i;
        ub |= unique_partners;
    }
    for (Mask t = ua; t; t &= t - 1) {
        if (repeated & (ub << trailing_zero(t))) return 0;
    }
    return 2;
}

int main(int argc, char** argv) {
    if (argc != 6) throw std::invalid_argument("m n cap outer_begin outer_end (0 = full end)");
    const int m = std::stoi(argv[1]), n = std::stoi(argv[2]), cap = std::stoi(argv[3]);
    const auto begin = std::stoull(argv[4]), end_arg = std::stoull(argv[5]);
    if (m < 1 || n < m || m > 126 || n > 126 || m + n > 126 || cap < 2 || cap > 9) throw std::invalid_argument("unsupported range");
    const auto as = supports(m, cap), bs = supports(n, cap);
    const auto end = end_arg ? end_arg : as.size();
    if (begin >= end || end > as.size()) throw std::invalid_argument("outer range");
    std::array<std::uint64_t, 3> counts{};
    std::vector<std::array<Mask, 2>> unknown;
    for (auto i = begin; i < end; ++i) {
        const auto& a = as[i];
        for (const auto& b : bs) {
            if (!canonical(a, b, m == n)) continue;
            const int status = classify(a.points, b.points, m, n);
            ++counts[status];
            if (status == 2) unknown.push_back({a.points, b.points});
        }
    }
    std::cout << "{\"m\":" << m << ",\"n\":" << n << ",\"cap\":" << cap
              << ",\"begin\":" << begin << ",\"end\":" << end << ",\"outer\":" << as.size()
              << ",\"inner\":" << bs.size() << ",\"raw\":" << (end-begin)*bs.size()
              << ",\"canonical\":" << counts[0]+counts[1]+counts[2]
              << ",\"excluded\":" << counts[0] << ",\"fair\":" << counts[1]
              << ",\"unknown\":" << counts[2] << ",\"unknown_masks\":[";
    for (std::size_t i = 0; i < unknown.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '[' << decimal(unknown[i][0]) << ',' << decimal(unknown[i][1]) << ']';
    }
    std::cout << "]}\n";
}
