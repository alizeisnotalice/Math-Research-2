// Exact directed-support peak inclusion-exclusion. U192 limb operations
// are inherited from r045-construction/cut_region_rootwise.cpp.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using Mask = std::uint64_t;
using U128 = __uint128_t;

struct U192 {
    std::array<std::uint64_t, 3> limbs{};

    U192& operator+=(const U192& other) {
        U128 carry = 0;
        for (int i = 0; i < 3; ++i) {
            const U128 total = U128(limbs[i]) + other.limbs[i] + carry;
            limbs[i] = std::uint64_t(total);
            carry = total >> 64;
        }
        assert(carry == 0);
        return *this;
    }

    U192& operator*=(unsigned scale) {
        U128 carry = 0;
        for (auto& limb : limbs) {
            const U128 total = U128(limb) * scale + carry;
            limb = std::uint64_t(total);
            carry = total >> 64;
        }
        assert(carry == 0);
        return *this;
    }

    unsigned divide(unsigned divisor) {
        U128 remainder = 0;
        for (int i = 2; i >= 0; --i) {
            const U128 total = (remainder << 64) | limbs[i];
            limbs[i] = std::uint64_t(total / divisor);
            remainder = total % divisor;
        }
        return unsigned(remainder);
    }

    U192& operator-=(const U192& other) {
        U128 borrow = 0;
        for (int i = 0; i < 3; ++i) {
            const U128 subtrahend = U128(other.limbs[i]) + borrow;
            const U128 minuend = limbs[i];
            limbs[i] = std::uint64_t(minuend - subtrahend);
            borrow = minuend < subtrahend;
        }
        assert(borrow == 0);
        return *this;
    }

    std::string decimal() const {
        U192 quotient = *this;
        std::string result;
        do {
            result.push_back(char('0' + quotient.divide(10)));
        } while (quotient.limbs[0] || quotient.limbs[1] || quotient.limbs[2]);
        std::reverse(result.begin(), result.end());
        return result;
    }
};
int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const int n = std::atoi(argv[1]), mid = std::atoi(argv[2]);
    const std::string widths = argv[3];
    if (n < 4 || n > 7 || mid < 2 || mid >= n-1 || int(widths.size()) != n) return 2;
    std::vector<int> vertices;
    for (int y = 0; y < n; ++y) {
        const int width = widths[y] - '0';
        if ((y%2 && width != mid) || (!(y%2) && (width < mid-1 || width > mid+1))) return 2;
        for (int x = 0; x < width; ++x) vertices.push_back(y*n+x);
    }
    const int count = int(vertices.size());
    if (count > 36) return 2;
    std::vector<Mask> support(count);
    const int band = (n+3)/4;
    const auto core = [&](int v) {
        return v/n%2 || v%n < mid-1;
    };
    const auto safe = [&](int v) {
        // The root is safe for every possible assignment of the flexible pairs.
        if (v%n >= band || !core(v)) return false;
        for (int u = 0; u < n*n; ++u) {
            if (std::abs(u/n-v/n)+std::abs(u%n-v%n) != 1) continue;
            if (u%n >= mid || !core(u)) return false;
        }
        return true;
    };
    Mask allowed = 0;
    for (int i = 0; i < count; ++i) {
        const int v = vertices[i];
        if (safe(v)) allowed |= Mask(1)<<i;
        for (int j = 0; j < count; ++j) {
            const int u = vertices[j];
            if (std::abs(u/n-v/n)+std::abs(u%n-v%n) != 1) continue;
            if ((core(v) || v%n == mid-1) && core(u)) support[i] |= Mask(1)<<j;
            if (!core(v) && v%n == mid && u == v-1) support[i] |= Mask(1)<<j;
        }
    }

    std::vector<Mask> mutual(count);
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < count; ++j)
            if ((support[i] & (Mask(1)<<j)) && (support[j] & (Mask(1)<<i)))
                mutual[i] |= Mask(1)<<j;
    std::vector<std::vector<Mask>> sets(count+1);
    std::size_t total = 0;
    const auto enumerate = [&](const auto& self, Mask candidates, Mask chosen) -> void {
        if (!candidates) {
            if (++total > 16'000'000) std::exit(3);
            sets[std::popcount(chosen)].push_back(chosen);
            return;
        }
        const int i = std::countr_zero(candidates);
        const Mask bit = Mask(1)<<i;
        candidates ^= bit;
        self(self,candidates,chosen);
        self(self,candidates & ~mutual[i],chosen|bit);
    };
    enumerate(enumerate,(Mask(1)<<count)-1,0);
    std::unordered_map<Mask,std::uint32_t> index;
    index.reserve(total);
    index.emplace(0,0);
    std::vector<U192> moments(total);
    moments[0].limbs[0] = 1;
    for (int j=2;j<=count;++j) moments[0] *= unsigned(j);
    std::uint32_t next = 1;
    U192 positive, negative;
    for (int k=1;k<=count;++k) {
        for (Mask s: sets[k]) {
            Mask neighbors = 0;
            for (Mask t=s;t;t&=t-1) neighbors |= support[std::countr_zero(t)];
            U192 moment;
            for (Mask t=s & ~neighbors;t;t&=t-1)
                moment += moments[index.at(s^(Mask(1)<<std::countr_zero(t)))];
            const auto remainder = moment.divide(std::popcount(s | neighbors));
            assert(remainder == 0);
            index.emplace(s,next);
            moments[next++] = moment;
            moment *= std::popcount(s & allowed);
            if (k%2) positive += moment;
            else negative += moment;
        }
        std::vector<Mask>().swap(sets[k]);
    }
    assert(next == total);
    positive -= negative;
    std::cout << positive.decimal() << '\n';
}
