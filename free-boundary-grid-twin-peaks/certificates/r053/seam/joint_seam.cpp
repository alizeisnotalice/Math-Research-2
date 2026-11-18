// Exact joint peak moments for a protected rectangle and its marked right edge.
// Dense cardinality layers adapted from the ranked-layer round-49 counter.
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
constexpr int MAX_AREA = 42;
constexpr int MAX_N = 6;
constexpr int MAX_PROTECTED = 36;
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
    Region(int width,int height) : w(width), n(height-1) {
        if (w<2 || w>6 || height<2 || height>7) throw std::invalid_argument("unsupported bottom rectangle");
        area=w*height;
        protected_count=(w-1)*n;
        mask_count=ONE<<n;
        for (int v=0;v<area;++v) {
            const int x=v%w,y=v/w;
            const Word b=ONE<<v;
            if (y && x<w-1) protected_vertices|=b;
            else unprotected_vertices|=b;
            closed[v]=b;
            for (auto [dx,dy]:{std::pair{-1,0},{1,0},{0,-1},{0,1}}) {
                const int xx=x+dx,yy=y+dy;
                if (0<=xx && xx<w && 0<=yy && yy<height) closed[v]|=ONE<<(yy*w+xx);
            }
        }
        row_choices.resize(n+1);
        row_choices[0].push_back({0,0,{}});
        for (int y = 0; y < n; ++y) {
            const Word capacity = ONE << w;
            auto& choices = row_choices[y+1];
            std::vector<unsigned char> index(capacity);
            for (Word mask = 0; mask < capacity; ++mask) {
                if (mask & (mask << 1)) continue;
                index[mask] = static_cast<unsigned char>(choices.size());
                RowChoice c{mask, unsigned(std::popcount(mask)), {}};
                for (Word bits = mask; bits; bits &= bits-1) {
                    const unsigned x = std::countr_zero(bits);
                    c.deleted.push_back({static_cast<unsigned char>((y+1)*w+x), index[mask ^ (ONE << x)]});
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
                for (int k = 0; k <= w*n; ++k) {
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
        auto end = counts.begin()+w*n+1;
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

    int w, n;
    int area=0;
    int protected_count=0;
    Word protected_vertices=0, unprotected_vertices=0, mask_count=0;
    std::array<Word, MAX_AREA> closed{};
    std::vector<Word> histogram;

private:
    std::vector<std::vector<RowChoice>> row_choices;
    std::vector<std::vector<RowState>> rows;
};

struct Resources {
    Word max_bytes = 0;
    Word independent_sets = 0;
};

Resources resources(const Region& region) {
    Resources plan;
    Word prev_count = 0;
    for (Word count : region.histogram) {
        plan.independent_sets += count;
        plan.max_bytes = std::max(plan.max_bytes, 16*(prev_count+count));
        prev_count = count;
    }
    return plan;
}

struct LayerEntry { Word counts, unprotected_predecessors; };
static_assert(sizeof(LayerEntry)==16);

using Moment = std::array<Word,4>;
struct Result {
    std::vector<Moment> moments;
    std::vector<Word> processed;
    std::vector<double> seconds;
};

Result solve(const Region& region, Word prime, bool log, bool verify_ranks) {
    const Modular mod(prime);
    std::array<Word, MAX_AREA+1> inverse{};
    Word factorial = 1;
    for (int d=1;d<=region.area;++d) {
        inverse[d]=mod.inverse(d);
        if (mod.multiply(d,inverse[d])!=1) throw std::logic_error("invalid inverse");
        factorial=mod.multiply(factorial,d);
    }
    const Word inv_first=inverse[region.area];
    const Word inv_second=mod.multiply(inv_first,inverse[region.area-1]);
    Result result{std::vector<Moment>(region.histogram.size()*region.mask_count),{1},{0.0}};
    auto record = [&](unsigned k, Word mask, Word closed, Word value, Word p_prev, Word b_prev,
                      Word cross_prev, Word unsafe_unique) {
        const Word P=region.protected_vertices & ~closed;
        const Word B=region.unprotected_vertices & ~closed;
        const Word p=std::popcount(P),b=std::popcount(B);
        Word edges=0;
        for (Word left=B;left;left&=left-1) edges+=std::popcount(P & region.closed[std::countr_zero(left)]);
        const Word T=mod.add(mod.multiply(b,value),b_prev);
        const Word V=mod.add(mod.add(mod.multiply(p*b-edges,value),mod.multiply(b,p_prev)),
                             mod.add(mod.multiply(p,b_prev),cross_prev));
        const Word A=mod.add(mod.multiply(edges,value),unsafe_unique);
        auto& out=result.moments[k*region.mask_count+mask];
        out[0]=mod.add(out[0],value);
        out[1]=mod.add(out[1],mod.multiply(inv_first,T));
        out[2]=mod.add(out[2],mod.multiply(inv_second,V));
        out[3]=mod.add(out[3],mod.multiply(inv_second,A));
    };
    record(0,0,0,factorial,0,0,0,0);
    std::vector<LayerEntry> previous{{factorial,0}};
    for (unsigned size=1;size<region.histogram.size();++size) {
        const auto start=std::chrono::steady_clock::now();
        Word processed=0;
        std::vector<LayerEntry> current(region.histogram[size]);
        RowPath path{};
        auto visit=[&](const RowPath& choices) {
            std::array<Word,MAX_N+2> prefix_less{},suffix{};
            std::array<unsigned,MAX_N+1> remaining{};
            remaining[0]=size;
            unsigned prev=0;
            for (int y=0;y<region.n;++y) {
                const unsigned row=choices[y],total=remaining[y];
                suffix[y]=region.rank_term(y,prev,row,total);
                prefix_less[y+1]=prefix_less[y]+(total ? region.rank_term(y,prev,row,total-1):0);
                remaining[y+1]=total-region.choice(y,row).count;
                prev=row;
            }
            for (int y=region.n-1;y>=0;--y) suffix[y]+=suffix[y+1];
            if (verify_ranks && suffix[0]!=processed) throw std::logic_error("enumeration rank wrong");
            Word p_prev=0,b_prev=0,cross_prev=0,closed=0,two_cover=0,selected_protected=0,mark=0;
            std::array<Word,MAX_AREA> deletion_values;
            prev=0;
            for (int y=0;y<region.n;++y) {
                const unsigned row=choices[y];
                for (const Deleted& d:region.choice(y,row).deleted) {
                    const unsigned substitute=d.next_row_index;
                    Word subindex=prefix_less[y]+region.rank_term(y,prev,substitute,remaining[y]-1)+suffix[y+2];
                    if (y+1<region.n) subindex+=region.rank_term(y+1,substitute,choices[y+1],remaining[y+1]);
                    if (subindex>=previous.size()) throw std::logic_error("subindex outside preceding layer");
                    if (verify_ranks) {
                        RowPath after=choices;
                        after[y]=d.next_row_index;
                        if (subindex!=region.direct_rank(after,size-1)) throw std::logic_error("shortcut rank wrong");
                    }
                    const unsigned v=d.protected_index; // actual vertex in the full rectangle
                    const auto& x=previous[subindex];
                    const Word item=ONE<<v;
                    if (item & region.protected_vertices) {
                        deletion_values[v]=x.counts;
                        selected_protected|=item;
                        p_prev=mod.add(p_prev,x.counts);
                        cross_prev=mod.add(cross_prev,x.unprotected_predecessors);
                    } else {
                        mark|=ONE<<y;
                        b_prev=mod.add(b_prev,x.counts);
                    }
                    const Word neighbors=region.closed[v];
                    two_cover|=closed & neighbors;
                    closed|=neighbors;
                }
                prev=row;
            }
            const Word singly_covered=closed & ~two_cover & region.unprotected_vertices;
            Word weighted=0;
            for (Word b=selected_protected;b;b&=b-1) {
                const int v=std::countr_zero(b);
                const Word unique=std::popcount(region.closed[v] & singly_covered);
                if (unique) weighted=mod.add(weighted,mod.multiply(unique,deletion_values[v]));
            }
            const Word C=mod.multiply(mod.add(p_prev,b_prev),inverse[std::popcount(closed)]);
            current[processed++]={C,b_prev};
            record(size-std::popcount(mark),mark,closed,C,p_prev,b_prev,cross_prev,weighted);
        };
        region.enumerate(0,0,size,path,visit);
        if (processed!=region.histogram[size]) throw std::logic_error("wrong layer population");
        result.processed.push_back(processed);
        previous=std::move(current);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        result.seconds.push_back(seconds);
        if (log) std::cerr << "layer "<<size<<" states "<<processed<<" elapsed_seconds "<<seconds<<'\n';
    }
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
    if (argc < 4) {
        throw std::invalid_argument("usage: joint_seam WIDTH HEIGHT PRIME_INDEX [--plan] [--progress] [--verify-ranks] [--max-table-mib MIB]");
    }
    const Word width=unsigned_argument(argv[1]);
    const Word height=unsigned_argument(argv[2]);
    const Word index=unsigned_argument(argv[3]);
    if (width>6 || height>7) throw std::invalid_argument("bad rectangle dimensions");
    if (index >= PRIME.size()) throw std::invalid_argument("prime index must be 0..4");
    Word max_mib = 1024;
    bool plan_only = false;
    bool progress = false;
    bool verify_ranks = false;
    for (int i = 4; i < argc; ++i) {
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
    const Region region{int(width),int(height)};
    const Resources plan = resources(region);
    const bool fits = plan.max_bytes <= max_mib * MIB;
    if (!plan_only && !fits) {
        throw std::runtime_error("two-table budget refused: needs " + std::to_string(plan.max_bytes) + " bytes");
    }
    Result answer;
    if (!plan_only) answer = solve(region, PRIME[index], progress, verify_ranks);
    std::cout << std::setprecision(10)
              << "{\"schema\":\"a022-joint-protected-seam-moments-v1\",\"mode\":\"" << (plan_only ? "plan" : "count")
              << "\",\"width\":" << width << ",\"height\":" << height << ",\"area\":" << region.area
              << ",\"prime_index\":" << index << ",\"prime\":\"" << PRIME[index] << '\"'
              << ",\"protected_vertices\":" << region.protected_count << ",\"mask_count\":" << region.mask_count
              << ",\"independent_sets\":" << plan.independent_sets
              << ",\"max_table_bytes\":" << plan.max_bytes  << ",\"table_budget_mib\":" << max_mib
              << ",\"fits_budget\":" << (fits ? "true" : "false") << ",\"states_by_layer\":";
    array_json(region.histogram);
    std::cout << ",\"values_by_layer\":";
    array_json(region.histogram);
    if (!plan_only) {
        std::cout << ",\"processed_by_layer\":";
        array_json(answer.processed);
        std::cout << ",\"wall_seconds_by_layer\":";
        array_json(answer.seconds);
        std::cout << ",\"moments\":[";
        for (std::size_t i=0;i<answer.moments.size();++i) {
            if (i) std::cout << ',';
            const auto& x=answer.moments[i];
            std::cout<<'['<<x[0]<<','<<x[1]<<','<<x[2]<<','<<x[3]<<']';
        }
        std::cout<<']';
    }
    std::cout << "}\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
}
