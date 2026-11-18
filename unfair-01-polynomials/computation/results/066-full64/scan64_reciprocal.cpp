// Independent pseudo-remainder sieve in the negative reciprocal variable.
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int N = 64;
using Poly = std::array<int, N + 1>;

std::vector<int> free_exponents() {
    std::vector<int> powers;
    for (int k = 36; k < 60; ++k) {
        if (k != 38) powers.push_back(N - k);
    }
    if (powers.size() != 23) throw std::logic_error("free exponents");
    return powers;
}

template<int p>
int exact_degree(Poly a) {
    Poly b{};
    for (int i = 0; i <= N; ++i) {
        if (a[i] < 0) a[i] += p;
        if (i > 0) b[i - 1] = a[i] * i % p;
    }
    int da = N;
    int db = N - 1;
    while (db >= 0) {
        while (da >= db) {
            const int shift = da - db;
            const int former = a[da];
            const int next = b[db];
            for (int i = 0; i < shift; ++i) a[i] = next * a[i] % p;
            for (int i = shift; i < da; ++i) {
                a[i] = (next * a[i] - former * b[i - shift] + p * p) % p;
            }
            a[da] = 0;
            do { --da; } while (da >= 0 && a[da] == 0);
        }
        a.swap(b);
        const int temporary = da;
        da = db;
        db = temporary;
    }
    return da;
}

Poly prepare(std::uint64_t mask, const std::vector<int>& powers) {
    Poly a{};
    a[0] = 1;
    for (int i = 0; i < 36; ++i) {
        if (i != 7) a[N - i] = (i % 2) ? -1 : 1;
    }
    a[N - 38] = 1;
    for (std::size_t i = 0; i < powers.size(); ++i) {
        if ((mask >> i) & 1) a[powers[i]] = (powers[i] % 2) ? -1 : 1;
    }
    return a;
}

int main(int argc, char** argv) {
    const auto powers = free_exponents();
    const auto bound = std::uint64_t{1} << powers.size();
    if (argc == 3 && std::string(argv[1]) == "--inspect") {
        const auto mask = std::stoull(argv[2]);
        if (mask >= bound) throw std::invalid_argument("mask");
        const auto a = prepare(mask, powers);
        std::cout << "[" << exact_degree<163>(a) << "," << exact_degree<1279>(a) << "]\n";
        return 0;
    }
    if (argc != 3) throw std::invalid_argument("start stop");
    const auto start = std::stoull(argv[1]);
    const auto stop = std::stoull(argv[2]);
    if (start >= stop || stop > bound) throw std::invalid_argument("range");
    std::uint64_t previous = start ^ (start >> 1);
    auto g = prepare(previous, powers);
    std::uint64_t first_pass = 0;
    bool first = true;
    std::cout << "{\"start\":" << start << ",\"stop\":" << stop << ",\"hits\":[";
    for (std::uint64_t index = start; index < stop; ++index) {
        const std::uint64_t gray = index ^ (index >> 1);
        if (gray != previous) {
            const int bit = std::countr_zero(gray ^ previous);
            const int pos = powers[bit];
            g[pos] = (gray >> bit) & 1 ? ((pos % 2) ? -1 : 1) : 0;
        }
        previous = gray;
        const int e = exact_degree<1279>(g);
        if (e < 7) continue;
        ++first_pass;
        const int d = exact_degree<163>(g);
        if (d < 7) continue;
        if (!first) std::cout << ",";
        first = false;
        std::cout << "{\"mask\":" << index << ",\"gcd163\":" << d << ",\"gcd1279\":" << e << '}';
    }
    std::cout << "],\"first_prime_pass\":" << first_pass << "}\n";
    return 0;
}
