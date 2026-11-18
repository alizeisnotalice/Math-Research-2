// Exhaust all physical perfect matchings of F_2^4 and collapse identical direction multisets.
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <thread>
#include <vector>

constexpr unsigned POINTS = 16, BLOCKS = 8, SUBSPACES = 67, FULL_J = 255;
using Mask = uint16_t;
using U64 = uint64_t;
using Entry = std::pair<U64, U64>;

struct Pattern {
    U64 key{}, physical_count{}, witness{};
    bool verify_representative{};
};

struct Report {
    U64 patterns{}, matchings{}, QJ{}, proper_def{}, density{}, equality{}, full_rank{};
    U64 represented_QJ{}, represented_proper_def{}, represented_density{}, represented_eq{};
    U64 direct_GF2{}, direct_walsh{};
    std::array<U64, 5> ambient_rank_counts{};
};

static void enumerate(Mask remaining, U64 directions, U64 physical, unsigned depth, std::vector<Entry>& data) {
    if (!remaining) {
        data.emplace_back(directions, physical);
        return;
    }
    unsigned first = __builtin_ctz(remaining);
    Mask tail = remaining & (remaining - 1);
    for (Mask choices = tail; choices; choices &= choices - 1) {
        unsigned mate = __builtin_ctz(choices);
        unsigned dir = first ^ mate;
        U64 add = U64(1) << (4 * (dir - 1));
        U64 packed = U64(first | (mate << 4)) << (8 * depth);
        enumerate(tail ^ (1u << mate), directions + add, physical | packed, depth + 1, data);
    }
}

static Mask xor_translate(Mask points, unsigned v) {
    Mask image = 0;
    for (; points; points &= points - 1) image |= 1u << (__builtin_ctz(points) ^ v);
    return image;
}

static unsigned rank_add(unsigned v, std::array<unsigned, 13>& pivot) {
    while (v) {
        unsigned b = 31u - __builtin_clz(v);
        if (!pivot[b]) {
            pivot[b] = v;
            return 1;
        }
        v ^= pivot[b];
    }
    return 0;
}

struct Geometry {
    std::vector<Mask> masks;
    std::array<unsigned, 65536> id{};
    std::array<unsigned, SUBSPACES> dimension{}, size{};
    std::array<unsigned, POINTS> normal{};
    std::array<std::array<unsigned, SUBSPACES>, SUBSPACES> join{}, meet{};

    Geometry() {
        std::set<Mask> known{1};
        for (unsigned x = 1; x < POINTS; ++x) {
            auto before = known;
            for (Mask W : before) if (!(W & (1u << x))) known.insert(W | xor_translate(W, x));
        }
        masks.assign(known.begin(), known.end());
        assert(masks.size() == SUBSPACES);
        std::array<unsigned, 5> histogram{};
        for (unsigned i = 0; i < SUBSPACES; ++i) {
            id[masks[i]] = i;
            size[i] = __builtin_popcount(masks[i]);
            dimension[i] = __builtin_ctz(size[i]);
            assert(size[i] == 1u << dimension[i]);
            ++histogram[dimension[i]];
        }
        assert((histogram == std::array<unsigned, 5>{1, 15, 35, 15, 1}));
        for (unsigned a = 0; a < SUBSPACES; ++a) {
            for (unsigned b = 0; b < SUBSPACES; ++b) {
                Mask total = 0;
                for (Mask points = masks[a]; points; points &= points - 1) {
                    total |= xor_translate(masks[b], __builtin_ctz(points));
                }
                meet[a][b] = id[masks[a] & masks[b]];
                join[a][b] = id[total];
                assert(masks[join[a][b]] == total);
            }
        }
        for (unsigned q = 1; q < POINTS; ++q) {
            Mask points = 0;
            for (unsigned x = 0; x < POINTS; ++x) {
                if ((__builtin_popcount(x & q) & 1) == 0) points |= 1u << x;
            }
            normal[q] = id[points];
        }
    }
};

static void walsh_verify(const Pattern& pattern, const Geometry& geom, Report& report) {
    std::array<unsigned, POINTS> tag{};
    std::array<Mask, BLOCKS> actual_normals{};
    for (unsigned i = 0; i < BLOCKS; ++i) {
        unsigned first = (pattern.witness >> (8 * i)) & 15;
        unsigned mate = (pattern.witness >> (8 * i + 4)) & 15;
        tag[first] = tag[mate] = i;
        actual_normals[i] = geom.masks[geom.normal[first ^ mate]];
    }
    constexpr unsigned SIZE = POINTS << BLOCKS;
    std::array<int, SIZE> walsh{};
    for (unsigned x = 0; x < SIZE; ++x) walsh[x] = 1 - 2 * ((x >> (4 + tag[x & 15])) & 1);
    for (unsigned width = 1; width < SIZE; width <<= 1) {
        for (unsigned start = 0; start < SIZE; start += 2 * width) {
            for (unsigned x = start; x < start + width; ++x) {
                int a = walsh[x], b = walsh[x + width];
                walsh[x] = a + b;
                walsh[x + width] = a - b;
            }
        }
    }
    for (unsigned x = 0; x < SIZE; ++x) {
        unsigned high = x >> 4;
        bool present = high && !(high & (high - 1)) && ((actual_normals[__builtin_ctz(high)] >> (x & 15)) & 1);
        assert(present == bool(walsh[x]));
        if (present) assert(std::abs(walsh[x]) == SIZE / BLOCKS);
    }
    report.direct_walsh += SIZE;
}

static void test(const Pattern& p, const Geometry& g, Report& report) {
    std::array<unsigned, BLOCKS> normals{};
    unsigned row = 0, w = 0;
    for (unsigned d = 1; d < POINTS; ++d) {
        unsigned copies = (p.key >> (4 * (d - 1))) & 15;
        for (unsigned j = 0; j < copies; ++j) {
            normals[row++] = g.normal[d];
            w = g.join[w][g.normal[d]];
        }
    }
    assert(row == BLOCKS);
    unsigned R = g.dimension[w], delta = BLOCKS - R;
    assert(R == 3 || R == 4);
    ++report.patterns;
    report.matchings += p.physical_count;
    report.ambient_rank_counts[R] += p.physical_count;
    std::array<unsigned, FULL_J + 1> full_normal_sum{};
    for (unsigned J = 0; J < FULL_J; ++J) {
        if (J) full_normal_sum[J] = g.join[full_normal_sum[J & (J - 1)]][normals[__builtin_ctz(J)]];
        int deficiency = __builtin_popcount(J) - int(g.dimension[full_normal_sum[J]]);
        assert(deficiency <= int(delta) - 1);
        ++report.proper_def;
        report.represented_proper_def += p.physical_count;
    }
    if (p.verify_representative) walsh_verify(p, g, report);
    for (unsigned Q = 0; Q < SUBSPACES; ++Q) {
        if (g.meet[w][Q] != Q) continue;
        std::array<unsigned, BLOCKS> ui{};
        for (unsigned i = 0; i < BLOCKS; ++i) ui[i] = g.meet[Q][normals[i]];
        std::array<unsigned, FULL_J + 1> sums{}, sizes{};
        for (unsigned J = 0; J <= FULL_J; ++J) {
            if (J) {
                unsigned first = __builtin_ctz(J), prior = J & (J - 1);
                sums[J] = g.join[sums[prior]][ui[first]];
                sizes[J] = sizes[prior] + g.size[ui[first]];
            }
            unsigned j_size = __builtin_popcount(J);
            unsigned r = j_size + g.dimension[sums[J]];
            assert(r <= j_size + g.dimension[Q]);
            ++report.QJ;
            report.represented_QJ += p.physical_count;
            if (p.verify_representative) {
                std::array<unsigned, 13> basis{};
                unsigned direct = 0;
                for (unsigned tmp = J; tmp; tmp &= tmp - 1) {
                    unsigned i = __builtin_ctz(tmp);
                    for (Mask available = g.masks[ui[i]]; available; available &= available - 1) {
                        unsigned word = (1u << 12) | (1u << (4 + i)) | __builtin_ctz(available);
                        direct += rank_add(word, basis);
                    }
                }
                assert(r == direct);
                ++report.direct_GF2;
            }
            if (r == R + BLOCKS) {
                ++report.full_rank;
                assert(Q == w && J == FULL_J && sizes[J] == BLOCKS * BLOCKS);
                continue;
            }
            ++report.density;
            report.represented_density += p.physical_count;
            unsigned cap2 = BLOCKS * (r + delta - 1);
            assert(2 * sizes[J] <= cap2);
            if (2 * sizes[J] == cap2) {
                ++report.equality;
                report.represented_eq += p.physical_count;
            }
        }
    }
}

int main(int argc, char **argv) {
    unsigned threads = argc == 2 ? std::stoul(argv[1]) : 12;
    if (!threads) return 2;
    const auto start = std::chrono::steady_clock::now();
    std::vector<Entry> physical;
    physical.reserve(2027025);
    enumerate(0xffff, 0, 0, 0, physical);
    U64 expected = 1;
    for (unsigned i = 1; i < POINTS; i += 2) expected *= i;
    assert(physical.size() == expected);
    std::sort(physical.begin(), physical.end());
    std::vector<Pattern> patterns;
    std::set<std::vector<unsigned>> multiplicities;
    for (size_t i = 0; i < physical.size();) {
        size_t j = i + 1;
        while (j < physical.size() && physical[i].first == physical[j].first) ++j;
        std::vector<unsigned> multiplicity;
        for (unsigned d = 0; d < 15; ++d) {
            unsigned count = (physical[i].first >> (4 * d)) & 15;
            if (count) multiplicity.push_back(count);
        }
        std::sort(multiplicity.begin(), multiplicity.end());
        patterns.push_back({physical[i].first, j - i, physical[i].second, multiplicities.insert(multiplicity).second});
        i = j;
    }
    Geometry geom;
    std::vector<Report> outputs(threads);
    std::vector<std::thread> pool;
    std::atomic<size_t> next{0};
    for (unsigned i = 0; i < threads; ++i) {
        pool.emplace_back([&geom, &patterns, &outputs, &next, i] {
            while (true) {
                size_t j = next.fetch_add(1);
                if (j >= patterns.size()) return;
                test(patterns[j], geom, outputs[i]);
            }
        });
    }
    for (auto& worker : pool) worker.join();
    Report all;
    for (const auto& r : outputs) {
        all.patterns += r.patterns;
        all.matchings += r.matchings;
        all.QJ += r.QJ;
        all.proper_def += r.proper_def;
        all.density += r.density;
        all.equality += r.equality;
        all.full_rank += r.full_rank;
        all.represented_QJ += r.represented_QJ;
        all.represented_proper_def += r.represented_proper_def;
        all.represented_density += r.represented_density;
        all.represented_eq += r.represented_eq;
        all.direct_GF2 += r.direct_GF2;
        all.direct_walsh += r.direct_walsh;
        for (unsigned k = 0; k < 5; ++k) all.ambient_rank_counts[k] += r.ambient_rank_counts[k];
    }
    assert(all.matchings == expected && all.full_rank == patterns.size());
    auto finish = std::chrono::steady_clock::now();
    std::cout << "{\"physical_matchings\":" << all.matchings
              << ",\"distinct_direction_multisets\":" << all.patterns
              << ",\"normal_rank_3_matchings\":" << all.ambient_rank_counts[3]
              << ",\"normal_rank_4_matchings\":" << all.ambient_rank_counts[4]
              << ",\"Q_J_compressed_checks\":" << all.QJ
              << ",\"Q_J_represented_checks\":" << all.represented_QJ
              << ",\"proper_subfamily_compressed_deficiency_checks\":" << all.proper_def
              << ",\"proper_subfamily_represented_deficiency_checks\":" << all.represented_proper_def
              << ",\"nonfull_rank_compressed_checks\":" << all.density
              << ",\"nonfull_rank_represented_checks\":" << all.represented_density
              << ",\"nonfull_rank_represented_equalities\":" << all.represented_eq
              << ",\"direct_GF2_candidates\":" << all.direct_GF2
              << ",\"direct_walsh_coefficients\":" << all.direct_walsh
              << ",\"direct_representatives\":" << multiplicities.size()
              << ",\"workers\":" << threads << ",\"wall_seconds\":"
              << std::chrono::duration<double>(finish-start).count() << "}\n";
}
