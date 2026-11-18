// Exact protected-peak histogram with side-maximum status modulo a 61-bit prime.
// Compact region masks and on-demand cardinality layers for n<=11, |H|<=66.
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
static_assert(sizeof(std::size_t) >= 8);
constexpr Word ONE = 1;
constexpr int MAX_AREA = 66;
constexpr int MAX_N = 11;
constexpr int MAX_PROTECTED = 55;
constexpr Word MIB = ONE << 20;
constexpr std::array<Word, 5> PRIME = {
    2305843009213693951ULL, 2305843009213693921ULL, 2305843009213693907ULL, 2305843009213693723ULL, 2305843009213693693ULL,
};

class Modular {
public:
    explicit Modular(Word modulus) : p(modulus), correction((ONE << 61) - p) {
        if (correction == 0 || correction > 259) throw std::invalid_argument("unsupported modulus");
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
    unsigned char row_index;
    unsigned char remaining;
};

struct Deleted {
    unsigned char protected_index;
    unsigned char next_row_index;
};

struct RowChoice {
    Word bits;
    unsigned count;
    std::vector<Deleted> deleted;
};

struct RowState {
    std::array<Word, MAX_PROTECTED + 1> counts{};
    std::array<std::vector<Branch>, MAX_PROTECTED + 1> branches;
    std::array<std::vector<Word>, MAX_PROTECTED + 1> ranks;
};

using RowPath = std::array<unsigned char, MAX_N>;

class Region {
public:
    explicit Region(std::string_view profile) : n(int(profile.size())) {
        if (n < 3 || n > MAX_N) throw std::invalid_argument("profile length must be 3..11");
        std::vector<int> widths(n), positions(n+1), site(n*n, -1), seam(n*n, -1);
        for (int y = 0; y < n; ++y) {
            const char c = profile[y];
            const int w = c-'0';
            if (c < '1' || c > '9' || w >= n) throw std::invalid_argument("invalid row width");
            if (y && std::abs(w-widths[y-1]) > 1) throw std::invalid_argument("non-Lipschitz profile");
            widths[y] = w;
            area += w;
            positions[y+1] = positions[y]+w-1;
            for (int x = 0; x < w-1; ++x) site[y*n+x] = positions[y]+x;
            seam[y*n+w-1] = y;
        }
        protected_count = area-n;
        if (area > MAX_AREA || protected_count > MAX_PROTECTED) throw std::invalid_argument("region exceeds capacity");
        all_protected = (ONE << protected_count)-1;
        all_seam = (ONE << n)-1;
        for (int v = 0; v < n*n; ++v) {
            const int k = site[v];
            if (k < 0 && seam[v] < 0) continue;
            if (k >= 0) protected_closed[k] |= ONE << k;
            for (auto [dx,dy] : {std::pair{-1,0},{1,0},{0,-1},{0,1}}) {
                const int x = v%n+dx, y = v/n+dy;
                if (x < 0 || x >= n || y < 0 || y >= n) continue;
                const int p = site[y*n+x], s = seam[y*n+x];
                if (k >= 0) {
                    if (p >= 0) protected_closed[k] |= ONE << p;
                    else if (s >= 0) seam_closed[k] |= ONE << s;
                    else throw std::logic_error("protected vertex has a cross-interface neighbor");
                } else if (p >= 0) seam_neighbors[seam[v]] |= ONE << p;
            }
        }
        row_choices.resize(n+1);
        row_choices[0].push_back({0,0,{}});
        for (int y = 0; y < n; ++y) {
            const Word capacity = ONE << (widths[y]-1);
            auto& choices = row_choices[y+1];
            std::vector<unsigned char> index(capacity);
            for (Word mask = 0; mask < capacity; ++mask) {
                if (mask & (mask << 1)) continue;
                index[mask] = static_cast<unsigned char>(choices.size());
                RowChoice c{mask, unsigned(std::popcount(mask)), {}};
                for (Word bits = mask; bits; bits &= bits-1) {
                    const unsigned x = std::countr_zero(bits);
                    c.deleted.push_back({static_cast<unsigned char>(positions[y]+x), index[mask ^ (ONE << x)]});
                }
                choices.push_back(std::move(c));
            }
        }
        rows.resize(n+1);
        for (int y = n; y >= 0; --y) {
            rows[y].resize(row_choices[y].size());
            for (std::size_t prev = 0; prev < rows[y].size(); ++prev) {
                RowState& state = rows[y][prev];
                if (y == n) {
                    state.counts[0] = 1;
                    continue;
                }
                const auto& next_rows = row_choices[y+1];
                for (int k = 0; k <= protected_count; ++k) {
                    auto& rank = state.ranks[k];
                    rank.resize(next_rows.size());
                    Word preceding = 0;
                    for (unsigned next = 0; next < next_rows.size(); ++next) {
                        rank[next] = preceding;
                        const auto& option = next_rows[next];
                        if (option.bits & row_choices[y][prev].bits || option.count > unsigned(k)) continue;
                        const unsigned left = k-option.count;
                        const Word choices = rows[y+1][next].counts[left];
                        if (choices == 0) continue;
                        preceding += choices;
                        state.branches[k].push_back({static_cast<unsigned char>(next), static_cast<unsigned char>(left)});
                    }
                    state.counts[k] = preceding;
                }
            }
        }
        const auto& counts = rows[0][0].counts;
        auto end = counts.begin()+protected_count+1;
        while (end != counts.begin() && end[-1] == 0) --end;
        histogram.assign(counts.begin(),end);
        if (histogram.empty() || histogram[0] != 1) throw std::logic_error("missing empty set");
    }

    template<class Visitor>
    void enumerate(int y, unsigned prev, unsigned left, RowPath& path, Visitor& visit) const {
        if (y == n) {
            visit(path);
            return;
        }
        for (const Branch& c : rows[y][prev].branches[left]) {
            path[y] = c.row_index;
            enumerate(y+1,c.row_index,c.remaining,path,visit);
        }
    }

    const Word& rank_term(int y,unsigned prev,unsigned row,unsigned left) const {
        return rows[y][prev].ranks[left][row];
    }

    const RowChoice& choice(int y,unsigned row) const { return row_choices[y+1][row]; }

    Word direct_rank(const RowPath& path, unsigned rem) const {
        unsigned prev = 0;
        Word rank = 0;
        for (int y = 0; y < n; ++y) {
            const unsigned row = path[y];
            rank += rank_term(y,prev,row,rem);
            rem -= choice(y,row).count;
            prev = row;
        }
        if (rem != 0) throw std::logic_error("bad path size");
        return rank;
    }

    int n;
    int area = 0;
    int protected_count = 0;
    Word all_protected = 0;
    Word all_seam = 0;
    std::array<Word, MAX_PROTECTED> protected_closed{}, seam_closed{};
    std::array<Word, MAX_N> seam_neighbors{};
    std::vector<Word> histogram;

private:
    std::vector<std::vector<RowChoice>> row_choices;
    std::vector<std::vector<RowState>> rows;
};

struct Resources {
    Word max_bytes = 0;
    Word legacy_hash_bytes = 0;
    Word independent_sets = 0;
};

Resources resources(const Region& region) {
    Resources plan;
    Word prev_count = 0, prev_slots = 0;
    for (Word count : region.histogram) {
        const Word slots = count ? std::bit_ceil((4*count+2)/3) : 0;
        plan.independent_sets += count;
        plan.max_bytes = std::max(plan.max_bytes, 8*(prev_count+count));
        plan.legacy_hash_bytes = std::max(plan.legacy_hash_bytes,16*(prev_slots+slots));
        prev_count = count;
        prev_slots = slots;
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

Result solve(const Region& region, Word prime, bool log, bool verify_ranks) {
    const Modular mod(prime);
    std::array<Word, MAX_AREA+1> inverse{};
    Word factorial = 1;
    for (int d = 1; d <= region.area; ++d) {
        inverse[d] = mod.inverse(d);
        if (mod.multiply(d,inverse[d]) != 1) throw std::logic_error("bad modular inverse");
        factorial = mod.multiply(factorial,d);
    }
    const auto inv_area = inverse[region.area];
    const Word inv_distinct = mod.multiply(inv_area,inverse[region.area-1]);
    auto bad_second_weight = [&](Word P, Word B) {
        const Word safe_out = region.all_protected & ~P;
        const Word unsafe_out = region.all_seam & ~B;
        Word result = Word(region.protected_count+std::popcount(B))*std::popcount(unsafe_out);
        for (Word bits = unsafe_out; bits; bits &= bits-1) {
            result -= std::popcount(region.seam_neighbors[std::countr_zero(bits)] & safe_out);
        }
        return result;
    };
    std::vector<Word> moments{factorial};
    std::vector<Word> first{mod.multiply(region.n,mod.multiply(factorial,inv_area))};
    std::vector<Word> second{mod.multiply(bad_second_weight(0,0),mod.multiply(factorial,inv_distinct))};
    Result result{{},{},{},{1},{0.0}};
    std::vector<Word> previous{factorial};
    for (unsigned size = 1; size < region.histogram.size(); ++size) {
        const auto start = std::chrono::steady_clock::now();
        Word processed = 0, moment = 0, first_sum = 0, second_sum = 0;
        std::vector<Word> current(region.histogram[size]);
        RowPath path{};
        auto visit = [&](const RowPath& choices) {
            std::array<Word, MAX_N+2> prefix_less{}, suffix{};
            std::array<unsigned, MAX_N+1> remaining{};
            remaining[0] = size;
            unsigned prev = 0;
            for (int y = 0; y < region.n; ++y) {
                const unsigned row = choices[y], total = remaining[y];
                suffix[y] = region.rank_term(y,prev,row,total);
                prefix_less[y+1] = prefix_less[y] + (total ? region.rank_term(y,prev,row,total-1) : 0);
                remaining[y+1] = total-region.choice(y,row).count;
                prev = row;
            }
            for (int y = region.n-1; y >= 0; --y) suffix[y] += suffix[y+1];
            if (verify_ranks && suffix[0] != processed) throw std::logic_error("enumeration order disagrees with rank");
            Word total = 0, P = 0, B = 0;
            prev = 0;
            for (int y = 0; y < region.n; ++y) {
                const unsigned row = choices[y];
                for (const Deleted& d : region.choice(y,row).deleted) {
                    const unsigned substitute = d.next_row_index;
                    Word subindex = prefix_less[y]+region.rank_term(y,prev,substitute,remaining[y]-1)+suffix[y+2];
                    if (y+1 < region.n) subindex += region.rank_term(y+1,substitute,choices[y+1],remaining[y+1]);
                    if (subindex >= previous.size()) throw std::logic_error("subindex outside preceding layer");
                    if (verify_ranks) {
                        RowPath alternate = choices;
                        alternate[y] = d.next_row_index;
                        if (subindex != region.direct_rank(alternate,size-1)) throw std::logic_error("shortcut rank wrong");
                    }
                    total = mod.add(total,previous[subindex]);
                    P |= region.protected_closed[d.protected_index];
                    B |= region.seam_closed[d.protected_index];
                }
                prev = row;
            }
            const int count = std::popcount(P)+std::popcount(B);
            const Word value = mod.multiply(total,inverse[count]);
            current[processed++] = value;
            moment = mod.add(moment,value);
            const int outside = region.n-std::popcount(B);
            if (outside) {
                first_sum = mod.add(first_sum,mod.multiply(outside,value));
                second_sum = mod.add(second_sum,mod.multiply(bad_second_weight(P,B),value));
            }
        };
        region.enumerate(0,0,size,path,visit);
        if (processed != region.histogram[size]) throw std::logic_error("wrong cardinality count");
        moments.push_back(moment);
        first.push_back(mod.multiply(first_sum,inv_area));
        second.push_back(mod.multiply(second_sum,inv_distinct));
        result.processed.push_back(processed);
        previous = std::move(current);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        result.seconds.push_back(elapsed);
        if (log) std::cerr << "layer " << size << " states " << processed << " elapsed_seconds " << elapsed << '\n';
    }
    auto invert = [&](const std::vector<Word>& sums) {
        std::vector<Word> counts(sums.size());
        for (std::size_t k = 0; k < sums.size(); ++k) {
            Word coefficient = 1;
            for (std::size_t j = 0; j <= k; ++j) {
                const Word term = mod.multiply(coefficient,sums[k]);
                counts[j] = (k-j)%2 ? mod.subtract(counts[j],term) : mod.add(counts[j],term);
                coefficient = coefficient*(k-j)/(j+1);
            }
        }
        return counts;
    };
    result.peak_frequency_residues = invert(moments);
    result.unprotected_max_frequency_residues = invert(first);
    result.bad_second_frequency_residues = invert(second);
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
        throw std::invalid_argument("usage: protected_mod PROFILE PRIME_INDEX [--plan] [--progress] [--verify-ranks] [--max-table-mib MIB]");
    }
    const std::string_view profile(argv[1]);
    const Word index = unsigned_argument(argv[2]);
    if (index >= PRIME.size()) throw std::invalid_argument("prime index must be 0..4");
    Word max_mib = 1024;
    bool plan_only = false;
    bool progress = false;
    bool verify_ranks = false;
    for (int i = 3; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        if (arg == "--plan") plan_only = true;
        else if (arg == "--progress") progress = true;
        else if (arg == "--verify-ranks") verify_ranks = true;
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
    if (!plan_only) answer = solve(region, PRIME[index], progress, verify_ranks);
    std::cout << std::setprecision(10)
              << "{\"schema\":\"a022-dense-protected-peaks-mod-v1\",\"mode\":\"" << (plan_only ? "plan" : "count")
              << "\",\"profile\":\"" << profile << "\",\"n\":" << region.n << ",\"area\":" << region.area
              << ",\"prime_index\":" << index << ",\"prime\":\"" << PRIME[index] << '\"'
              << ",\"protected_vertices\":" << region.protected_count
              << ",\"independent_sets\":" << plan.independent_sets
              << ",\"max_table_bytes\":" << plan.max_bytes << ",\"legacy_open_address_payload_bytes\":" << plan.legacy_hash_bytes << ",\"table_budget_mib\":" << max_mib
              << ",\"fits_budget\":" << (fits ? "true" : "false") << ",\"states_by_layer\":";
    array_json(region.histogram);
    std::cout << ",\"values_by_layer\":";
    array_json(region.histogram);
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
