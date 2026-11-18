// Independent replay: combinations by mask successor, lexicographic vectors,
// and explicit counts of representations for each sum.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Word = std::uint64_t;
using Points = std::vector<int>;
struct Support {
    Word mask;
    Points points, reverse;
    int gcd;
};

bool precedes(const Points& a, const Points& b) {
    return a.size() == b.size() ? a < b : a.size() < b.size();
}

std::vector<Support> generate(int degree, int cap) {
    std::vector<Support> result;
    const Word limit = Word{1} << (degree - 1);
    for (int r = 0; r <= std::min(cap - 2, degree - 1); ++r) {
        std::vector<Support> size_group;
        for (Word x = (Word{1} << r) - 1; x < limit;) {
            Points p{0};
            int g = degree;
            for (Word z = x; z; z &= z - 1) {
                int i = std::countr_zero(z) + 1;
                p.push_back(i);
                g = std::gcd(g, i);
            }
            p.push_back(degree);
            Points q;
            for (auto it = p.rbegin(); it != p.rend(); ++it) q.push_back(degree - *it);
            size_group.push_back({(x << 1) | (Word{1} << degree) | 1, p, q, g});
            if (r == 0) break;
            const Word low = x & -x;
            const Word carry = x + low;
            x = carry | ((carry ^ x) >> (std::countr_zero(low) + 2));
        }
        std::sort(size_group.begin(), size_group.end(), [](const auto& u, const auto& v) {return u.points < v.points;});
        result.insert(result.end(), size_group.begin(), size_group.end());
    }
    return result;
}

bool keep(const Support& a, const Support& b, bool same) {
    if (std::gcd(a.gcd, b.gcd) != 1 || (same && precedes(b.points, a.points))) return false;
    const auto* left = &a.reverse;
    const auto* right = &b.reverse;
    if (same && precedes(*right, *left)) std::swap(left, right);
    if (a.points != *left) return precedes(a.points, *left);
    return !precedes(*right, b.points);
}

int decide(const Points& a, const Points& b) {
    std::array<unsigned char, 64> multiplicity{};
    bool collision = false;
    for (int i : a) for (int j : b) {
        collision |= ++multiplicity[i + j] > 1;
    }
    if (!collision) return 1;
    std::array<bool, 64> ua{}, ub{};
    ua[0] = ua[a.back()] = ub[0] = ub[b.back()] = true;
    for (int i : a) for (int j : b) {
        if (multiplicity[i + j] == 1) ua[i] = ub[j] = true;
    }
    for (int i : a) if (ua[i]) for (int j : b) {
        if (ub[j] && multiplicity[i + j] > 1) return 0;
    }
    return 2;
}

int main(int argc, char** argv) {
    if (argc != 6) throw std::invalid_argument("m n cap start end");
    const int m = std::stoi(argv[1]), n = std::stoi(argv[2]), cap = std::stoi(argv[3]);
    const auto lo = std::stoull(argv[4]), endpoint = std::stoull(argv[5]);
    if (m < 1 || n < m || m + n > 62 || cap < 2 || cap > 9) throw std::invalid_argument("unsupported range");
    const auto a = generate(m, cap), b = generate(n, cap);
    const auto hi = endpoint ? endpoint : a.size();
    if (hi > a.size() || lo >= hi) throw std::invalid_argument("outer range");
    std::array<Word, 3> ct{};
    std::vector<std::array<Word, 2>> uncertain;
    for (auto i = lo; i < hi; ++i) for (const auto& v : b) {
        if (!keep(a[i], v, m == n)) continue;
        const int s = decide(a[i].points, v.points);
        ++ct[s];
        if (s == 2) uncertain.push_back({a[i].mask, v.mask});
    }
    std::cout << "{\"m\":" << m << ",\"n\":" << n << ",\"cap\":" << cap
              << ",\"begin\":" << lo << ",\"end\":" << hi << ",\"outer\":" << a.size()
              << ",\"inner\":" << b.size() << ",\"raw\":" << (hi-lo)*b.size()
              << ",\"canonical\":" << ct[0]+ct[1]+ct[2] << ",\"excluded\":" << ct[0]
              << ",\"fair\":" << ct[1] << ",\"unknown\":" << ct[2] << ",\"unknown_masks\":[";
    for (std::size_t i = 0; i < uncertain.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '[' << uncertain[i][0] << ',' << uncertain[i][1] << ']';
    }
    std::cout << "]}\n";
}
