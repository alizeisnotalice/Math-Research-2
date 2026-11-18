// Exhaustively sieve the 23 free coefficients of a degree-64 binary product.
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

constexpr int degree = 64;
constexpr int threshold = 7;
using Polynomial = std::array<int, degree + 1>;
constexpr std::array<int, 23> positions = {
    36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59
};
constexpr std::uint64_t bound = std::uint64_t{1} << positions.size();

template<int p>
constexpr auto inverses() {
    std::array<int, p> values{};
    values[1] = 1;
    for (int i = 2; i < p; ++i) {
        values[i] = p - ((p / i) * values[p % i]) % p;
    }
    return values;
}

template<int p>
int gcd_degree(Polynomial a, bool shortcut) {
    constexpr auto inv = inverses<p>();
    Polynomial b{};
    for (int i = 1; i <= degree; ++i) {
        b[i - 1] = i * a[i] % p;
    }
    Polynomial* u = &a;
    Polynomial* v = &b;
    int n = degree;
    int m = degree - 1;
    while (m >= 0 && (!shortcut || m >= threshold)) {
        const int multiplier = inv[(*v)[m]];
        for (int top = n; top >= m; --top) {
            if (!(*u)[top]) continue;
            const int shift = top - m;
            const int factor = (*u)[top] * multiplier % p;
            for (int i = 0; i < m; ++i) {
                (*u)[shift + i] = ((*u)[shift + i] + p * p - factor * (*v)[i]) % p;
            }
            (*u)[top] = 0;
        }
        int remainder_degree = m - 1;
        while (remainder_degree >= 0 && !(*u)[remainder_degree]) {
            --remainder_degree;
        }
        auto* temp = u;
        u = v;
        v = temp;
        n = m;
        m = remainder_degree;
    }
    return m < 0 ? n : m;
}

Polynomial make_product(std::uint64_t mask) {
    Polynomial a{};
    for (int i = 0; i <= 35; ++i) a[i] = i != 7;
    a[38] = a[64] = 1;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        a[positions[i]] = (mask >> i) & 1;
    }
    return a;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--inspect") {
        const auto mask = std::stoull(argv[2]);
        if (mask >= bound) throw std::invalid_argument("mask");
        const auto a = make_product(mask);
        std::cout << "[" << gcd_degree<163>(a, false) << "," << gcd_degree<1279>(a, false) << "]\n";
        return 0;
    }
    if (argc != 3) throw std::invalid_argument("start stop");
    const auto start = std::stoull(argv[1]);
    const auto stop = std::stoull(argv[2]);
    if (start >= stop || stop > bound) throw std::invalid_argument("range");
    std::uint64_t first_pass = 0;
    std::cout << "{\"start\":" << start << ",\"stop\":" << stop << ",\"hits\":[";
    bool first = true;
    for (std::uint64_t mask = start; mask < stop; ++mask) {
        const auto a = make_product(mask);
        const int d = gcd_degree<163>(a, true);
        if (d < threshold) continue;
        ++first_pass;
        const int e = gcd_degree<1279>(a, true);
        if (e < threshold) continue;
        if (!first) std::cout << ",";
        first = false;
        std::cout << "{\"mask\":" << mask << ",\"gcd163\":" << d << ",\"gcd1279\":" << e << '}';
    }
    std::cout << "],\"first_prime_pass\":" << first_pass << "}\n";
    return 0;
}
