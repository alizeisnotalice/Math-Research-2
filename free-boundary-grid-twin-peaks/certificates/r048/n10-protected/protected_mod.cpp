// Exact protected-peak histogram with side-maximum status modulo a 61-bit prime.
// Compact region masks and on-demand cardinality layers for n<=10, |H|<=50.
#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using Word = std::uint64_t;
using Wide = __uint128_t;
constexpr Word ONE = 1;
constexpr int MAX_AREA = 50;
constexpr Word MIB = ONE << 20;
constexpr std::array<Word, 4> PRIME = {
    2305843009213693951ULL, 2305843009213693921ULL, 2305843009213693907ULL, 2305843009213693723ULL,
};

class Modular {
public:
    explicit Modular(Word modulus) : p(modulus), correction((ONE << 61) - p) {
        if (correction == 0 || correction > 229) throw std::invalid_argument("unsupported modulus");
    }

    Word add(Word a, Word b) const {
        const Word sum = a + b;
        return sum >= p ? sum - p : sum;
    }

    Word subtract(Word a, Word b) const {
        return a >= b ? a - b : a + p - b;
    }

    Word multiply(Word a, Word b) const {
        constexpr Word low = (ONE << 61) - 1;
        const Wide product = Wide(a) * b;
        // Since 2^61 = correction mod p, two folds leave less than 2p.
        const Wide first = (Word(product) & low) + (product >> 61) * correction;
        const Word second = (Word(first) & low) + Word(first >> 61) * correction;
        return second >= p ? second - p : second;
    }

    Word inverse(Word a) const {
        Word power = p - 2;
        Word result = 1;
        while (power) {
            if (power & 1) result = multiply(result, a);
            a = multiply(a, a);
            power >>= 1;
        }
        return result;
    }

private:
    Word p;
    Word correction;
};

struct Branch {
    Word placed;
    unsigned char row_index;
    unsigned char remaining;
};

struct RowState {
    std::array<Word, MAX_AREA + 1> counts{};
    std::array<std::vector<Branch>, MAX_AREA + 1> branches;
};

class Region {
public:
    explicit Region(std::string_view profile) : n(int(profile.size())) {
        if (n < 3 || n > 10) throw std::invalid_argument("profile length must be 3..10");
        std::vector<int> widths;
        widths.reserve(n);
        for (char digit : profile) {
            const int width = digit - '0';
            if (width < 1 || width >= n) throw std::invalid_argument("invalid row width");
            widths.push_back(width);
            area += width;
        }
        if (area > MAX_AREA) throw std::invalid_argument("region exceeds 50 vertices");

        for (int y = 1; y < n; ++y) {
            if (std::abs(widths[y] - widths[y - 1]) > 1) {
                throw std::invalid_argument("row widths must be 1-Lipschitz");
            }
        }
        std::vector<int> offset(n + 1), local(n * n, -1);
        std::vector<std::vector<Word>> valid_rows(n + 1, {0});
        for (int y = 0; y < n; ++y) {
            offset[y + 1] = offset[y] + widths[y];
            for (int x = 0; x < widths[y]; ++x) {
                local[y * n + x] = offset[y] + x;
                global.push_back(y * n + x);
            }
            for (Word row = 1; row < (ONE << (widths[y] - 1)); ++row) {
                if (!(row & (row << 1))) valid_rows[y + 1].push_back(row);
            }
        }

        for (int pos = 0; pos < area; ++pos) {
            const int v = global[pos];
            closed[pos] = ONE << pos;
            bool safe = true;
            for (auto [dx, dy] : {std::pair{-1, 0}, {1, 0}, {0, -1}, {0, 1}}) {
                const int x = v % n + dx;
                const int y = v / n + dy;
                if (x < 0 || x >= n || y < 0 || y >= n) continue;
                const int u = local[y * n + x];
                if (u < 0) safe = false;
                else closed[pos] |= ONE << u;
            }
            if (safe != (v % n < widths[v / n] - 1)) throw std::logic_error("protected geometry mismatch");
            if (!safe) unprotected |= ONE << pos;
        }

        rows.resize(n + 1);
        for (int y = n; y >= 0; --y) {
            rows[y].resize(valid_rows[y].size());
            for (std::size_t prev = 0; prev < valid_rows[y].size(); ++prev) {
                RowState& state = rows[y][prev];
                if (y == n) {
                    state.counts[0] = 1;
                    continue;
                }
                for (std::size_t next = 0; next < valid_rows[y + 1].size(); ++next) {
                    const Word selected = valid_rows[y + 1][next];
                    if (selected & valid_rows[y][prev]) continue;
                    const int chosen = std::popcount(selected);
                    const Word placed = selected << offset[y];
                    for (int left = 0; left <= area - chosen; ++left) {
                        const Word suffix_count = rows[y + 1][next].counts[left];
                        if (suffix_count == 0) continue;
                        state.counts[left + chosen] += suffix_count;
                        state.branches[left + chosen].push_back({placed, static_cast<unsigned char>(next),
                                                                 static_cast<unsigned char>(left)});
                    }
                }
            }
        }

        const auto& counts = rows.front().front().counts;
        auto end = counts.begin() + area + 1;
        while (end != counts.begin() && end[-1] == 0) --end;
        histogram.assign(counts.begin(), end);
        if (histogram.empty() || histogram[0] != 1) throw std::logic_error("missing empty set");
    }

    template<class Visitor>
    void enumerate(int y, unsigned prev, unsigned remaining, Word mask, Visitor& visit) const {
        if (y == n) {
            visit(mask);
            return;
        }
        for (const Branch& branch : rows[y][prev].branches[remaining]) {
            enumerate(y + 1, branch.row_index, branch.remaining, mask | branch.placed, visit);
        }
    }

    int n;
    int area = 0;
    Word unprotected = 0;
    std::array<Word, MAX_AREA> closed{};
    std::vector<int> global;
    std::vector<Word> histogram;

private:
    std::vector<std::vector<RowState>> rows;
};

std::size_t slots_for(Word entries) {
    if (entries == 0) return 0;
    // Power-of-two capacity with occupancy at most 3/4.
    return static_cast<std::size_t>(std::bit_ceil((4 * entries + 2) / 3));
}

struct Entry { Word key = 0; Word residue = 0; };
static_assert(sizeof(Entry) == 16 && sizeof(std::size_t) >= 8);

class Layer {
public:
    explicit Layer(std::size_t capacity = 0) : table(capacity), mask(capacity ? capacity - 1 : 0) {}

    Word operator[](Word key) const {
        if (key == 0) return empty_set;
        std::size_t at = hash(key) & mask;
        while (table[at].key != key) {
            if (table[at].key == 0) throw std::logic_error("submoment missing");
            at = (at + 1) & mask;
        }
        return table[at].residue;
    }

    void insert(Word key, Word value) {
        std::size_t at = hash(key) & mask;
        while (table[at].key) {
            if (table[at].key == key) throw std::logic_error("duplicate independent set");
            at = (at + 1) & mask;
        }
        table[at] = {key, value};
    }

    Word empty_set = 0;  // The only cardinality-zero key is held outside the table.

private:
    static Word hash(Word key) {
        key = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
        key = (key ^ (key >> 27)) * 0x94d049bb133111ebULL;
        return key ^ (key >> 31);
    }

    std::vector<Entry> table;
    std::size_t mask;
};

struct Resources {
    std::vector<std::size_t> slots;
    std::size_t max_bytes = 0;
    Word independent_sets = 0;
};

Resources resources(const Region& region) {
    Resources plan;
    for (Word count : region.histogram) {
        plan.slots.push_back(slots_for(count));
        plan.independent_sets += count;
    }
    // The empty-set residue sits in a scalar, not in a hash slot.
    plan.slots[0] = 0;
    for (std::size_t k = 1; k < plan.slots.size(); ++k) {
        plan.max_bytes = std::max(plan.max_bytes, sizeof(Entry) * (plan.slots[k - 1] + plan.slots[k]));
    }
    return plan;
}

struct Result {
    std::vector<Word> peak_frequency_residues;
    std::vector<Word> unprotected_max_frequency_residues;
    std::vector<Word> bad_second_frequency_residues;
    std::vector<Word> processed;
    std::vector<double> seconds;
};

Result solve(const Region& region, const Resources& plan, Word prime, bool progress) {
    const Modular mod(prime);
    std::array<Word, MAX_AREA + 1> inverse{};
    Layer previous;
    previous.empty_set = 1;
    for (int d = 1; d <= region.area; ++d) {
        inverse[d] = mod.inverse(Word(d));
        if (mod.multiply(d, inverse[d]) != 1) throw std::logic_error("invalid inverse");
        previous.empty_set = mod.multiply(previous.empty_set, Word(d));
    }
    const Word full_factorial = previous.empty_set;
    const Word prev_factorial = mod.multiply(full_factorial, inverse[region.area]);
    std::vector<Word> moments{full_factorial};
    std::vector<Word> max_moments{mod.multiply(std::popcount(region.unprotected), prev_factorial)};
    const Word protected_vertices = ((ONE << region.area) - 1) & ~region.unprotected;
    const Word ordered_top_pair_inverse = mod.multiply(inverse[region.area], inverse[region.area - 1]);
    auto bad_second_weight = [&](Word neighborhood) {
        const Word unsafe_outside = region.unprotected & ~neighborhood;
        const Word safe_outside = protected_vertices & ~neighborhood;
        Word weight = Word(std::popcount(safe_outside) + std::popcount(neighborhood)) * std::popcount(unsafe_outside);
        for (Word bits = unsafe_outside; bits; bits &= bits - 1) {
            weight -= std::popcount(region.closed[std::countr_zero(bits)] & safe_outside);
        }
        return weight;
    };
    std::vector<Word> second_moments{mod.multiply(mod.multiply(ordered_top_pair_inverse, full_factorial),
                                                bad_second_weight(0))};
    Result result{{}, {}, {}, {1}, {0.0}};
    for (std::size_t size = 1; size < region.histogram.size(); ++size) {
        const auto start = std::chrono::steady_clock::now();
        Layer current(plan.slots[size]);
        Word processed = 0;
        Word moment_sum = 0;
        Word unprotected_max_moment_sum = 0;
        Word bad_second_moment_sum = 0;
        auto visit = [&](Word set) {
            Word total = 0;
            Word neighborhood = 0;
            for (Word remaining = set; remaining; remaining &= remaining - 1) {
                const unsigned position = std::countr_zero(remaining);
                total = mod.add(total, previous[set ^ (ONE << position)]);
                neighborhood |= region.closed[position];
            }
            const Word value = mod.multiply(total, inverse[std::popcount(neighborhood)]);
            current.insert(set, value);
            moment_sum = mod.add(moment_sum, value);
            const auto outside = std::popcount(region.unprotected & ~neighborhood);
            if (outside) {
                unprotected_max_moment_sum = mod.add(unprotected_max_moment_sum, mod.multiply(outside, value));
                bad_second_moment_sum = mod.add(bad_second_moment_sum, mod.multiply(bad_second_weight(neighborhood), value));
            }
            ++processed;
        };
        region.enumerate(0, 0, static_cast<unsigned>(size), 0, visit);
        if (processed != region.histogram[size]) throw std::logic_error("wrong layer population");
        result.processed.push_back(processed);
        moments.push_back(moment_sum);
        max_moments.push_back(mod.multiply(inverse[region.area], unprotected_max_moment_sum));
        second_moments.push_back(mod.multiply(ordered_top_pair_inverse, bad_second_moment_sum));
        previous = std::move(current);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        result.seconds.push_back(elapsed);
        if (progress) {
            std::cerr << "layer " << size << " states " << processed << " elapsed_seconds " << elapsed << '\n';
        }
    }
    auto invert = [&](const std::vector<Word>& sums) {
        std::vector<Word> counts(sums.size(), 0);
        for (std::size_t k = 0; k < sums.size(); ++k) {
            Word coefficient = 1;
            for (std::size_t j = 0; j <= k; ++j) {
                const Word term = mod.multiply(coefficient, sums[k]);
                counts[j] = (k - j) % 2 ? mod.subtract(counts[j], term) : mod.add(counts[j], term);
                coefficient = coefficient * (k - j) / (j + 1);
            }
        }
        return counts;
    };
    result.peak_frequency_residues = invert(moments);
    result.unprotected_max_frequency_residues = invert(max_moments);
    result.bad_second_frequency_residues = invert(second_moments);
    return result;
}

Word unsigned_argument(std::string_view value) {
    Word number = 0;
    const auto [ptr, error] = std::from_chars(value.data(), value.data() + value.size(), number);
    if (error != std::errc{} || ptr != value.data() + value.size()) {
        throw std::invalid_argument("expected an unsigned integer argument");
    }
    return number;
}

template<class T>
void array_json(const std::vector<T>& values) {
    std::cout << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << values[i];
    }
    std::cout << ']';
}

int main(int argc, char** argv) try {
    if (argc < 3) {
        throw std::invalid_argument("usage: protected_mod PROFILE PRIME_INDEX [--plan] [--progress] [--max-table-mib MIB]");
    }
    const std::string_view profile(argv[1]);
    const Word index = unsigned_argument(argv[2]);
    if (index >= PRIME.size()) throw std::invalid_argument("prime index must be 0..3");
    Word max_mib = 1024;
    bool plan_only = false;
    bool progress = false;
    for (int i = 3; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        if (arg == "--plan") plan_only = true;
        else if (arg == "--progress") progress = true;
        else if (arg == "--max-table-mib" && i + 1 < argc) max_mib = unsigned_argument(argv[++i]);
        else throw std::invalid_argument("invalid option");
    }
    if (max_mib == 0 || max_mib > std::numeric_limits<Word>::max() / MIB) {
        throw std::invalid_argument("table budget out of range");
    }
    const Region region(profile);
    const Resources plan = resources(region);
    const bool fits = plan.max_bytes <= max_mib * MIB;
    if (!plan_only && !fits) {
        throw std::runtime_error("two-table budget refused: needs " + std::to_string(plan.max_bytes) + " bytes");
    }
    Result answer;
    if (!plan_only) answer = solve(region, plan, PRIME[index], progress);
    std::cout << std::setprecision(10)
              << "{\"schema\":\"a022-protected-peaks-mod-v2\",\"mode\":\"" << (plan_only ? "plan" : "count")
              << "\",\"profile\":\"" << profile << "\",\"n\":" << region.n << ",\"area\":" << region.area
              << ",\"prime_index\":" << index << ",\"prime\":\"" << PRIME[index] << '\"'
              << ",\"protected_vertices\":" << region.area - std::popcount(region.unprotected)
              << ",\"independent_sets\":" << plan.independent_sets
              << ",\"max_table_bytes\":" << plan.max_bytes << ",\"table_budget_mib\":" << max_mib
              << ",\"fits_budget\":" << (fits ? "true" : "false") << ",\"states_by_layer\":";
    array_json(region.histogram);
    std::cout << ",\"slots_by_layer\":";
    array_json(plan.slots);
    if (!plan_only) {
        std::cout << ",\"processed_by_layer\":";
        array_json(answer.processed);
        std::cout << ",\"wall_seconds_by_layer\":";
        array_json(answer.seconds);
        std::cout << ",\"peak_frequency_residues\":";
        array_json(answer.peak_frequency_residues);
        std::cout << ",\"unprotected_max_frequency_residues\":";
        array_json(answer.unprotected_max_frequency_residues);
        std::cout << ",\"bad_second_frequency_residues\":";
        array_json(answer.bad_second_frequency_residues);
    }
    std::cout << "}\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
