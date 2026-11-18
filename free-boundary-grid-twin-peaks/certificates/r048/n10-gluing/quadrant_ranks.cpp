// Exact connector ranks in connected-growth orders on P_m square P_m, m<=5.
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

struct RootRanks {
    std::vector<std::vector<Wide>> by_connector;
    std::uint64_t positive_prefixes;
};

RootRanks ranks_for(int root, const Grid& grid, const Suffix& suffix) {
    const Mask root_bit = Mask(1) << root;
    const Mask below_root = root_bit - 1;
    const Mask compressed_full = grid.full >> 1;
    std::vector<Count> prefix(compressed_full + 1);
    prefix[0] = 1;  // The root bit is implicit in every state.
    std::array<Mask, 25> compressed_bit{};
    for (int v = 0; v < grid.area; ++v) {
        if (v != root) compressed_bit[v] = Mask(1) << (v - int(v > root));
    }
    RootRanks result{std::vector<std::vector<Wide>>(grid.m, std::vector<Wide>(grid.area+1)), 0};
    if (root / grid.m == grid.m-1) {
        result.by_connector[root % grid.m][1] = suffix.orders[root_bit];
    }
    for (Mask packed = 0; packed <= compressed_full; ++packed) {
        const Count value = prefix[packed];
        if (!value) continue;
        ++result.positive_prefixes;
        const Mask S = (packed & below_root) | root_bit | ((packed & ~below_root) << 1);
        const int next_rank = 1 + std::popcount(S);
        for (Mask available = grid.frontier(S); available; available &= available-1) {
            const int v = std::countr_zero(available);
            Count& next_prefix = prefix[packed | compressed_bit[v]];
            if (MAX_COUNT - next_prefix < value) throw std::overflow_error("prefix count exceeds uint64");
            next_prefix += value;
            if (v / grid.m != grid.m-1) continue;
            const Count extension = suffix.orders[S | (Mask(1) << v)];
            if (extension == MAX_COUNT) throw std::overflow_error("connected extension count was saturated");
            result.by_connector[v % grid.m][next_rank] += Wide(value)*extension;
        }
    }
    if (prefix[compressed_full] != suffix.orders[root_bit]) throw std::logic_error("root order total mismatch");
    for (const auto& ranks : result.by_connector) {
        Wide sum = 0;
        for (Wide count : ranks) sum += count;
        if (sum != prefix[compressed_full]) throw std::logic_error("rank histogram does not sum to rooted total");
    }
    return result;
}

int main(int argc, char** argv) try {
    if (argc != 2 || std::string(argv[1]).size() != 1 || argv[1][0] < '2' || argv[1][0] > '5') {
        throw std::invalid_argument("usage: quadrant_ranks SIDE, where SIDE is 2, 3, 4 or 5");
    }
    const int m = argv[1][0]-'0';
    const Grid grid(m);
    const Suffix suffix = suffix_counts(grid);
    std::cout << "{\"m\":" << m << ",\"saturated_suffix_masks\":" << suffix.saturated
              << ",\"largest_dense_payload_bytes\":"
              << sizeof(Count)*(std::uint64_t(grid.full+1)+std::uint64_t((grid.full>>1)+1))
              << ",\"rooted_counts\":[";
    for (int v = 0; v < grid.area; ++v) {
        if (v) std::cout << ',';
        const Count value = suffix.orders[Mask(1) << v];
        if (value == MAX_COUNT) throw std::overflow_error("single-root exact count exceeds uint64");
        std::cout << '"' << value << '"';
    }
    std::cout << "],\"roots\":[";
    bool first = true;
    for (int y = 0; y < m; ++y) {
        for (int x = 0; x < m-1 && 2*x <= m-1; ++x) {
            const int root = y*m+x;
            const RootRanks ranks = ranks_for(root, grid, suffix);
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"root\":" << root << ",\"positive_prefixes\":" << ranks.positive_prefixes << ",\"rank_counts\":[";
            for (int c = 0; c < m; ++c) {
                if (c) std::cout << ',';
                std::cout << '[';
                for (int j = 1; j <= grid.area; ++j) {
                    if (j > 1) std::cout << ',';
                    std::cout << '"' << decimal(ranks.by_connector[c][j]) << '"';
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
