// Independent replay: combinations by mask successor, lexicographic vectors,
// and explicit counts of representations for each sum.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <span>
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

using Word = __uint128_t;
using Points = std::span<const std::uint8_t>;
struct Support {
    std::array<std::uint8_t,9> points{}, reverse{};
    std::uint8_t terms=0, gcd=0;

    Points forward() const { return {points.data(),terms}; }
    Word bit_mask() const {
        Word out=0;
        for (int i : forward()) out |= Word{1}<<i;
        return out;
    }
};

bool precedes(int a_count, const std::array<std::uint8_t,9>& a,
              int b_count, const std::array<std::uint8_t,9>& b) {
    return a_count == b_count ? a < b : a_count < b_count;
}

std::size_t support_count(int degree, int cap) {
    std::size_t size=1, terms=1;
    for (int r=1; r<=std::min(cap-2,degree-1); ++r) {
        terms=terms*(degree-r)/r;
        size+=terms;
    }
    return size;
}

std::vector<Support> generate(int degree, int cap) {
    std::vector<Support> result;
    result.reserve(support_count(degree,cap));
    const Word limit = Word{1} << (degree - 1);
    for (int r=0; r<=std::min(cap-2,degree-1); ++r) {
        const auto from=result.size();
        for (Word x = (Word{1}<<r)-1; x<limit;) {
            Support p;
            p.terms=static_cast<std::uint8_t>(r+2);
            int i=1;
            int common=degree;
            for (Word z=x; z; z&=z-1) {
                const int j=trailing_zero(z)+1;
                p.points[i++]=static_cast<std::uint8_t>(j);
                common=std::gcd(common,j);
            }
            p.points[i]=static_cast<std::uint8_t>(degree);
            p.gcd=static_cast<std::uint8_t>(common);
            for (int j=0;j<p.terms;++j) p.reverse[j]=static_cast<std::uint8_t>(degree-p.points[p.terms-1-j]);
            result.push_back(p);
            if (!r) break;
            const Word low=x&-x,carry=x+low;
            x=carry|((carry^x)>>(trailing_zero(low)+2));
        }
        std::sort(result.begin()+from,result.end(),[](const Support& a,const Support& b){return a.points < b.points;});
    }
    return result;
}

bool keep(const Support& a,const Support& b,bool equal) {
    if (std::gcd(a.gcd,b.gcd)!=1 || (equal && precedes(b.terms,b.points,a.terms,a.points))) return false;
    auto* first=&a;
    auto* second=&b;
    if (equal && precedes(b.terms,b.reverse,a.terms,a.reverse)) std::swap(first,second);
    if (a.terms==first->terms && a.points==first->reverse) {
        return !precedes(second->terms,second->reverse,b.terms,b.points);
    }
    return precedes(a.terms,a.points,first->terms,first->reverse);
}

int decide(Points a, Points b) {
    std::array<unsigned char, 128> multiplicity{};
    bool collision = false;
    for (int i : a) for (int j : b) {
        collision |= ++multiplicity[i + j] > 1;
    }
    if (!collision) return 1;
    std::array<bool, 128> ua{}, ub{};
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
    if (m < 1 || n < m || m > 126 || n > 126 || m + n > 126 || cap < 2 || cap > 9) throw std::invalid_argument("unsupported range");
    const auto a = generate(m, cap), b = generate(n, cap);
    const auto hi = endpoint ? endpoint : a.size();
    if (hi > a.size() || lo >= hi) throw std::invalid_argument("outer range");
    std::array<Word, 3> ct{};
    std::vector<std::array<Word, 2>> uncertain;
    for (auto i = lo; i < hi; ++i) for (const auto& v : b) {
        if (!keep(a[i], v, m == n)) continue;
        const int s = decide(a[i].forward(), v.forward());
        ++ct[s];
        if (s == 2) uncertain.push_back({a[i].bit_mask(), v.bit_mask()});
    }
    std::cout << "{\"m\":" << m << ",\"n\":" << n << ",\"cap\":" << cap
              << ",\"begin\":" << lo << ",\"end\":" << hi << ",\"outer\":" << a.size()
              << ",\"inner\":" << b.size() << ",\"raw\":" << (hi-lo)*b.size()
              << ",\"canonical\":" << decimal(ct[0]+ct[1]+ct[2]) << ",\"excluded\":" << decimal(ct[0])
              << ",\"fair\":" << decimal(ct[1]) << ",\"unknown\":" << decimal(ct[2]) << ",\"unknown_masks\":[";
    for (std::size_t i = 0; i < uncertain.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '[' << decimal(uncertain[i][0]) << ',' << decimal(uncertain[i][1]) << ']';
    }
    std::cout << "]}\n";
}
