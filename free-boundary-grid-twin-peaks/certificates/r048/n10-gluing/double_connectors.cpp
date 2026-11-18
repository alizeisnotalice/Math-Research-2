// Joint last-connector ranks; grid and suffix DP adapted from quadrant_ranks.cpp.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using Mask = std::uint32_t;
using Count = std::uint64_t;
using Wide = __uint128_t;
constexpr Count MAX_COUNT = std::numeric_limits<Count>::max();

std::string decimal(Wide value) {
    std::string result;
    do {
        result.push_back(char('0' + value % 10));
        value /= 10;
    } while (value);
    std::reverse(result.begin(), result.end());
    return result;
}

class Grid {
public:
    explicit Grid(int side) : m(side), area(m*m), full((Mask(1) << area)-1), low_n(area/2),
                             low_mask((Mask(1) << low_n)-1), neighbors(area) {
        for (int v = 0; v < area; ++v) {
            int x = v % m, y = v / m;
            if (x) neighbors[v] |= Mask(1) << (v-1);
            if (x < m-1) neighbors[v] |= Mask(1) << (v+1);
            if (y) neighbors[v] |= Mask(1) << (v-m);
            if (y < m-1) neighbors[v] |= Mask(1) << (v+m);
        }
        low = table(0, low_n);
        high = table(low_n, area-low_n);
    }
    Mask frontier(Mask S) const {
        return (low[S & low_mask] | high[S >> low_n]) & ~S;
    }
    int m, area;
    Mask full;
private:
    std::vector<Mask> table(int base, int width) const {
        std::vector<Mask> result(Mask(1) << width);
        for (Mask s = 1; s < result.size(); ++s) {
            result[s] = result[s & (s-1)] | neighbors[base+std::countr_zero(s)];
        }
        return result;
    }
    int low_n;
    Mask low_mask;
    std::vector<Mask> neighbors, low, high;
};

struct Suffix {
    std::vector<Count> orders;
    std::uint64_t saturated = 0;
};

Suffix suffix_counts(const Grid& grid) {
    Suffix suffix{std::vector<Count>(grid.full+1), 0};
    suffix.orders[grid.full] = 1;
    for (Mask S = grid.full-1; S; --S) {
        Wide value = 0;
        for (Mask available = grid.frontier(S); available; available &= available-1) {
            value += suffix.orders[S | (available & -available)];
        }
        if (value >= MAX_COUNT) {
            ++suffix.saturated;
            suffix.orders[S] = MAX_COUNT;
        } else {
            suffix.orders[S] = Count(value);
        }
    }
    return suffix;
}

using Histogram = std::vector<std::vector<std::vector<Wide>>>;

Histogram double_ranks(int root, const Grid& grid, const Suffix& suffix) {
    const Mask root_bit = Mask(1) << root;
    const Mask below_root = root_bit - 1;
    const Mask compressed_full = grid.full >> 1;
    std::vector<Count> prefix(compressed_full+1);
    prefix[0] = 1;
    std::array<Mask, 25> compressed_bit{};
    Mask bottom_arc = 0;
    Mask right_arc = 0;
    for (int v = 0; v < grid.area; ++v) {
        if (v != root) compressed_bit[v] = Mask(1) << (v-int(v>root));
        if (v/grid.m == grid.m-1 && v%grid.m < grid.m-1) bottom_arc |= Mask(1) << v;
        if (v%grid.m == grid.m-1 && v/grid.m < grid.m-1) right_arc |= Mask(1) << v;
    }
    Histogram table(grid.m-1, std::vector<std::vector<Wide>>(grid.m-1, std::vector<Wide>(grid.area+1)));
    for (Mask packed = 0; packed <= compressed_full; ++packed) {
        const Count value = prefix[packed];
        if (!value) continue;
        const Mask S = (packed & below_root) | root_bit | ((packed & ~below_root) << 1);
        const int next_rank = 1 + std::popcount(S);
        for (Mask available = grid.frontier(S); available; available &= available-1) {
            const int v = std::countr_zero(available);
            const Mask bit = Mask(1) << v;
            Count& next_prefix = prefix[packed | compressed_bit[v]];
            if (MAX_COUNT-next_prefix < value) throw std::overflow_error("prefix count exceeds uint64");
            next_prefix += value;
            if (!(bit & (bottom_arc | right_arc))) continue;
            const Count extension = suffix.orders[S | bit];
            if (extension == MAX_COUNT) throw std::overflow_error("connected suffix saturated");
            const Wide words = Wide(value)*extension;
            if (bit & bottom_arc) {
                for (Mask other = S & right_arc; other; other &= other-1) {
                    table[v%grid.m][std::countr_zero(other)/grid.m][next_rank] += words;
                }
            } else {
                for (Mask other = S & bottom_arc; other; other &= other-1) {
                    table[std::countr_zero(other)%grid.m][v/grid.m][next_rank] += words;
                }
            }
        }
    }
    if (prefix[compressed_full] != suffix.orders[root_bit]) throw std::logic_error("root order total mismatch");
    for (const auto& row : table) {
        for (const auto& bins : row) {
            Wide sum = 0;
            for (Wide item : bins) sum += item;
            if (sum != suffix.orders[root_bit]) throw std::logic_error("double rank marginal mismatch");
        }
    }
    return table;
}

int main(int argc, char** argv) try {
    if (argc != 2 || std::string(argv[1]).size() != 1 || argv[1][0] < '2' || argv[1][0] > '5') {
        throw std::invalid_argument("usage: double_connectors SIDE (2..5)");
    }
    const Grid grid(argv[1][0]-'0');
    const Suffix suffix = suffix_counts(grid);
    std::cout << "{\"m\":" << grid.m << ",\"saturated_suffix_masks\":" << suffix.saturated << ",\"roots\":[";
    bool first = true;
    for (int y = 0; y < grid.m-1; ++y) {
        for (int x = 0; x < grid.m-1; ++x) {
            const int r = y*grid.m+x;
            const Histogram table = double_ranks(r, grid, suffix);
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"root\":" << r << ",\"last_connector_rank_counts\":[";
            for (int a = 0; a < grid.m-1; ++a) {
                if (a) std::cout << ',';
                std::cout << '[';
                for (int b = 0; b < grid.m-1; ++b) {
                    if (b) std::cout << ',';
                    std::cout << '[';
                    for (int j = 1; j <= grid.area; ++j) {
                        if (j > 1) std::cout << ',';
                        std::cout << '"' << decimal(table[a][b][j]) << '"';
                    }
                    std::cout << ']';
                }
                std::cout << ']';
            }
            std::cout << "]}";
        }
    }
    std::cout << "]}\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
