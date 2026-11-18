// Exact first-merge rank, low-branch size and next-contact delay, n <= 4.
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

using Mask = std::uint32_t;
using Count = std::uint64_t;
using Tail = std::array<Count, 17>;
using Row = std::array<Count, 16>;

class Counter {
  public:
    explicit Counter(int n)
        : n_(n), N_(n * n), full_((Mask{1} << N_) - 1),
          around_(full_ + 1), prefix_(full_ + 1), suffix_(full_ + 1) {
        for (int u = 0; u < N_; ++u) {
            for (int v = 0; v < N_; ++v) {
                if (distance(u, v) == 1) {
                    adjacency_[u] |= Mask{1} << v;
                }
            }
        }
        for (Mask set = 1; set <= full_; ++set) {
            const int u = first(set);
            const Mask rest = set ^ (Mask{1} << u);
            around_[set] = adjacency_[u] | around_[rest];
            if (!rest) {
                prefix_[set][u] = 1;
                continue;
            }
            for (Mask choices = set; choices; drop(choices)) {
                const int v = first(choices);
                const Mask smaller = set ^ (Mask{1} << v);
                if (adjacency_[v] & smaller) {
                    for (int root = 0; root < N_; ++root) {
                        prefix_[set][root] += prefix_[smaller][root];
                    }
                }
            }
        }
        suffix_[full_] = 1;
        for (Mask s = full_ - 1; s > 0; --s) {
            for (Mask frontier = around_[s] & ~s; frontier; drop(frontier)) {
                suffix_[s] += suffix_[s | (Mask{1} << first(frontier))];
            }
        }
    }

    void run(int shard, int shards) {
        std::map<std::tuple<int, int, int, int>, Count> sums;
        for (Mask s = static_cast<Mask>(shard); s < full_; s += shards) {
            if (std::popcount(s) < 2) {
                continue;
            }
            const Mask c = component(s, first(s));
            const Mask d = s ^ c;
            if (!d || component(d, first(d)) != d) {
                continue;
            }
            const Mask join = around_[c] & around_[d] & ~s;
            accumulate(c, d, join, sums);
            accumulate(d, c, join, sums);
        }
        for (const auto& [key, value] : sums) {
            const auto [i, b, h, distance] = key;
            std::cout << i << ' ' << b << ' ' << h << ' '
                      << distance << ' ' << value << '\n';
        }
    }

  private:
    static int first(Mask x) { return std::countr_zero(x); }
    static void drop(Mask& x) { x &= x - 1; }

    int distance(int u, int v) const {
        return std::abs(u % n_ - v % n_) + std::abs(u / n_ - v / n_);
    }

    Mask component(Mask s, int start) const {
        Mask seen = Mask{1} << start;
        for (Mask expanded = seen | (around_[seen] & s); expanded != seen;
             expanded = seen | (around_[seen] & s)) {
            seen = expanded;
        }
        return seen;
    }

    Tail until_contact(Mask s, Mask target) {
        assert((s & target) == 0);
        if (!target) {
            Tail result{};
            result[N_ + 1 - std::popcount(s)] = suffix_[s];
            return result;
        }
        auto key = (std::uint64_t(target) << N_) | s;
        if (const auto it = memo_.find(key); it != memo_.end()) {
            return it->second;
        }
        Tail result{};
        for (Mask frontier = around_[s] & ~s; frontier; drop(frontier)) {
            const Mask bit = Mask{1} << first(frontier);
            if (target & bit) {
                result[1] += suffix_[s | bit];
            } else {
                const auto later = until_contact(s | bit, target);
                for (int h = 1; h < N_; ++h) {
                    result[h + 1] += later[h];
                }
            }
        }
        memo_.emplace(key, result);
        return result;
    }

    static Count choose(int n, int k) {
        Count product = 1;
        for (int j = 1; j <= k; ++j) {
            product = product * (n - j + 1) / j;
        }
        return product;
    }

    void accumulate(Mask high, Mask low, Mask join,
                    std::map<std::tuple<int, int, int, int>, Count>& sums) {
        if (!join) {
            return;
        }
        const int b = std::popcount(low);
        const Mask s = high | low;
        const int i = std::popcount(s) + 1;
        const Count shuffle = choose(i - 2, b);
        std::array<Count, 8> pairs{};
        for (Mask high_f = high; high_f; drop(high_f)) {
            const int p = first(high_f);
            for (Mask low_f = low; low_f; drop(low_f)) {
                const int q = first(low_f);
                pairs[distance(p, q)] += prefix_[high][p] * prefix_[low][q];
            }
        }
        for (Mask w_f = join; w_f; drop(w_f)) {
            const Mask w = Mask{1} << first(w_f);
            const Mask next = s | w;
            const Mask target = around_[low] & ~next;
            const auto tail = until_contact(next, target);
            Count observed = 0;
            for (int h = 1; h <= N_; ++h) {
                observed += tail[h];
                if (tail[h]) {
                    for (int d = 2; d <= 2 * n_ - 2; ++d) {
                        if (pairs[d]) {
                            sums[{i, b, h, d}] += shuffle * pairs[d] * tail[h];
                        }
                    }
                }
            }
            assert(observed == suffix_[next]);
        }
    }

    int n_;
    int N_;
    Mask full_;
    std::array<Mask, 16> adjacency_{};
    std::vector<Mask> around_;
    std::vector<Row> prefix_;
    std::vector<Count> suffix_;
    std::unordered_map<std::uint64_t, Tail> memo_;
};

int main(int argc, char** argv) {
    assert(argc == 4);
    const int n = std::atoi(argv[1]);
    const int shard = std::atoi(argv[2]);
    const int shards = std::atoi(argv[3]);
    assert(n >= 2 && n <= 4 && 0 <= shard && shard < shards);
    Counter(n).run(shard, shards);
}
