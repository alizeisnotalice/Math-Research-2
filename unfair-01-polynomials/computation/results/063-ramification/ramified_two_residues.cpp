// Meet in the middle on the four linear equations for a binary product.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int prime1 = 163;
constexpr int prime2 = 1279;
constexpr int half = 22;
constexpr uint32_t range = uint32_t{1} << half;
using Payload = std::array<int, 4>;

int positive_mod(int value, int modulus) {
    int reduced = value % modulus;
    return reduced < 0 ? reduced + modulus : reduced;
}

uint64_t encode(const Payload& part) {
    const uint64_t offset = static_cast<uint64_t>(positive_mod(part[0], prime1));
    return ((offset * prime1 + positive_mod(part[1], prime1)) * prime2
            + positive_mod(part[2], prime2)) * prime2 + positive_mod(part[3], prime2);
}

std::array<Payload, 57> columns(int root1, int root) {
    std::array<Payload, 57> result{};
    int previous1 = 0;
    int power1 = 1;
    int previous = 0;
    int power = 1;
    for (int k = 0; k <= 56; ++k) {
        result[k] = {power1, (k*previous1) % prime1, power, (k*previous) % prime2};
        previous1 = power1;
        power1 = (power1*root1) % prime1;
        previous = power;
        power = (power*root) % prime2;
    }
    return result;
}

void accumulate(Payload& total, const Payload& term, int sign) {
    for (std::size_t i = 0; i < total.size(); ++i) {
        total[i] += sign*term[i];
    }
}

struct Record {
    uint64_t key;
    uint32_t mask;
};
}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        throw std::invalid_argument("pass residue roots modulo 163 and 1279");
    }
    const int root1 = positive_mod(std::stoi(argv[1]), prime1);
    const int root = std::stoi(argv[2]);
    if (root < 2 || root >= prime2) {
        throw std::invalid_argument("invalid residue root");
    }
    const auto coeff = columns(root1, root);
    Payload fixed{};
    for (int k = 0; k < 12; ++k) {
        if (k != 7) {
            accumulate(fixed, coeff[k], 1);
        }
    }
    accumulate(fixed, coeff[56], 1);

    std::vector<Record> left;
    left.reserve(range);
    Payload total{};
    uint32_t previous = 0;
    for (uint32_t i = 0; i < range; ++i) {
        const uint32_t code = i ^ (i >> 1);
        if (i != 0) {
            const uint32_t change = code ^ previous;
            accumulate(total, coeff[12+std::countr_zero(change)], (code & change) ? 1 : -1);
        }
        left.push_back({encode(total),code});
        previous = code;
    }
    std::sort(left.begin(), left.end(), [](const Record& a, const Record& b) {
        return a.key < b.key || (a.key == b.key && a.mask < b.mask);
    });

    total = fixed;
    for (int& term : total) {
        term = -term;
    }
    previous = 0;
    for (uint32_t i = 0; i < range; ++i) {
        const uint32_t code = i ^ (i >> 1);
        if (i != 0) {
            const uint32_t change = code ^ previous;
            accumulate(total, coeff[34+std::countr_zero(change)], (code & change) ? -1 : 1);
        }
        const uint64_t key = encode(total);
        auto row = std::lower_bound(left.begin(), left.end(), key, [](const Record& a, uint64_t b) {
            return a.key < b;
        });
        for (; row != left.end() && row->key == key; ++row) {
            std::cout << (static_cast<uint64_t>(code) << half) + row->mask << '\n';
        }
        previous = code;
    }
}
