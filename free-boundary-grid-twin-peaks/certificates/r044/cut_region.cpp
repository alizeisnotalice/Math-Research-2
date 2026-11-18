// Exact rooted connected-order counts. The independent-set recurrence
// follows a022-curved-single-region-n7-r076/source/curved_single_region_root.cpp.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
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

void enumerate(int n, const std::vector<std::vector<Mask>>& row_masks,
               std::vector<std::vector<Mask>>& by_size, int row = 0,
               Mask previous = 0, Mask full = 0) {
    if (row == n) {
        by_size[std::popcount(full)].push_back(full);
        return;
    }
    for (Mask candidate : row_masks[row]) {
        if ((candidate & previous) == 0) {
            enumerate(n, row_masks, by_size, row + 1, candidate,
                      full | (candidate << (n * row)));
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const int n = std::atoi(argv[1]);
    const std::string profile(argv[2]);
    if (n < 3 || n > 8 || int(profile.size()) != n) return 2;
    Mask region = 0;
    std::vector<std::vector<Mask>> row_masks(n);
    for (int y = 0; y < n; ++y) {
        const int width = profile[y] - '0';
        if (width < 1 || width >= n) return 2;
        region |= ((Mask(1) << width) - 1) << (n * y);
        for (Mask s = 0; s < (Mask(1) << width); ++s) {
            if ((s & (s << 1)) == 0) row_masks[y].push_back(s);
        }
    }
    const int vertices = std::popcount(region);
    // Both parity sums are bounded by 2^36 * 36 * 36! < 2^182.
    if (vertices > 36) return 2;
    std::vector<Mask> closed_neighborhood(n * n);
    Mask allowed = 0;
    for (int v = 0; v < n * n; ++v) {
        Mask closed = Mask(1) << v;
        for (int u = 0; u < n * n; ++u) {
            if (std::abs(u % n - v % n) + std::abs(u / n - v / n) == 1) {
                closed |= Mask(1) << u;
            }
        }
        closed_neighborhood[v] = closed & region;
        if (v % n < (n + 3) / 4 && !(closed & ~region)) allowed |= Mask(1) << v;
    }

    std::vector<std::vector<Mask>> by_size(vertices + 1);
    enumerate(n, row_masks, by_size);
    std::size_t total = 0;
    int alpha = 0;
    for (int k = 0; k <= vertices; ++k) {
        total += by_size[k].size();
        if (!by_size[k].empty()) alpha = k;
    }
    if (total > std::numeric_limits<std::uint32_t>::max()) return 2;
    std::unordered_map<Mask, std::uint32_t> index;
    index.reserve(total);
    index.emplace(0, 0);
    std::vector<U192> moment(total);
    moment[0].limbs[0] = 1;
    for (int i = 2; i <= vertices; ++i) moment[0] *= unsigned(i);
    std::uint32_t offset = 1;
    U192 positive, negative;
    for (int k = 1; k <= alpha; ++k) {
        for (Mask s : by_size[k]) {
            U192 sum;
            Mask closed = 0;
            for (Mask rem = s; rem; rem &= rem - 1) {
                const int v = std::countr_zero(rem);
                sum += moment[index.at(s ^ (Mask(1) << v))];
                closed |= closed_neighborhood[v];
            }
            const unsigned remainder = sum.divide(unsigned(std::popcount(closed)));
            assert(remainder == 0);
            index.emplace(s, offset);
            moment[offset++] = sum;
            sum *= unsigned(std::popcount(s & allowed));
            if (k % 2 == 1) positive += sum;
            else negative += sum;
        }
        std::vector<Mask>().swap(by_size[k]);
    }
    assert(offset == total && index.size() == total);
    positive -= negative;
    std::cout << "{\"n\":" << n << ",\"profile\":\"" << profile
              << "\",\"vertices\":" << vertices << ",\"allowed_vertices\":"
              << std::popcount(allowed) << ",\"alpha\":" << alpha
              << ",\"independent_sets\":" << total << ",\"rooted_count\":\""
              << positive.decimal() << "\"}\n";
}
