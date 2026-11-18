// Counts one-core-peak orders with directed support for flexible vertices.
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U = __uint128_t;
std::string decimal(U value) {
    std::string s;
    do { s.push_back(char('0' + value % 10)); value /= 10; } while (value);
    std::reverse(s.begin(), s.end());
    return s;
}

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
    if (count > 25) return 2;
    std::vector<std::uint32_t> support(count);
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
    const std::uint32_t cap = 1u << count;
    std::vector<U> dp(cap);
    for (int i = 0; i < count; ++i) {
        const int v = vertices[i];
        if (safe(v)) dp[1u<<i] = 1;
        for (int j = 0; j < count; ++j) {
            const int u = vertices[j];
            if (std::abs(u/n-v/n)+std::abs(u%n-v%n) != 1) continue;
            if ((core(v) || v%n == mid-1) && core(u)) support[i] |= 1u<<j;
            if (!core(v) && v%n == mid && u == v-1) support[i] |= 1u<<j;
        }
    }
    for (std::uint32_t s = 1; s < cap; ++s) {
        if (std::has_single_bit(s)) continue;
        U total = 0;
        for (std::uint32_t t = s; t; t &= t-1) {
            const unsigned i = std::countr_zero(t);
            const std::uint32_t prev = s ^ (1u<<i);
            if (prev & support[i]) total += dp[prev];
        }
        dp[s] = total;
    }
    std::cout << decimal(dp.back()) << '\n';
}
