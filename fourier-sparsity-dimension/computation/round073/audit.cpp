#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <thread>
#include <unordered_map>
#include <vector>

using Word = std::uint64_t;
using Count = std::uint64_t;

int parity(int x) { return std::popcount(static_cast<unsigned>(x)) & 1; }
int low(Word x) { return std::countr_zero(x); }

Word translate(Word input, int offset) {
    Word output = 0;
    for (; input; input &= input - 1) output |= Word(1) << (low(input) ^ offset);
    return output;
}

Word append(Word basis, int vector) { return basis | translate(basis, vector); }

Word span(Word points) {
    Word space = 1;
    for (; points; points &= points - 1) {
        int next = low(points);
        if (!(space & (Word(1) << next))) space = append(space, next);
    }
    return space;
}

int dimension(Word subspace) { return std::countr_zero(static_cast<unsigned>(std::popcount(subspace))); }

struct Split {
    Word space;
    Word kernel;
    int dim;
    std::vector<int> address;
    std::vector<int> physical_kernel;
    std::vector<int> bands;
    std::vector<Word> frequency_masks;
    std::unordered_map<Word, std::vector<Word>> activities;
};

struct Geometry {
    int n;
    int size;
    std::vector<Split> splits;

    explicit Geometry(int n) : n(n), size(1 << n) {
        std::vector<std::vector<Word>> by_rank = {{1}};
        for (int dim = 1; dim <= n; ++dim) {
            std::vector<Word> next;
            for (Word W : by_rank.back())
                for (int v = 1; v < size; ++v)
                    if (!(W & (Word(1) << v))) next.push_back(append(W, v));
            std::sort(next.begin(), next.end());
            next.erase(std::unique(next.begin(), next.end()), next.end());
            by_rank.push_back(std::move(next));
        }
        for (const auto& subspaces : by_rank) {
            for (Word H : subspaces) {
                Split split{};
                split.space = H;
                split.dim = dimension(H);
                for (int t = 0; t < size; ++t) {
                    bool orthogonal = true;
                    for (Word points = H; points; points &= points - 1)
                        orthogonal &= !parity(t & low(points));
                    if (orthogonal) {
                        split.kernel |= Word(1) << t;
                        split.physical_kernel.push_back(t);
                    }
                    const Word band = translate(H, t);
                    if (t == low(band)) {
                        split.bands.push_back(t);
                        split.frequency_masks.push_back(band);
                    }
                }
                for (int t = 0; t < size; ++t) {
                    if (t == low(translate(split.kernel, t))) split.address.push_back(t);
                }
                assert(split.address.size() == static_cast<unsigned>(std::popcount(H)));
                for (int dim = 0; dim <= split.dim; ++dim) {
                    for (Word W : by_rank[dim]) {
                        if ((W & H) != W) continue;
                        std::map<Word, Word> classes;
                        for (unsigned i = 0; i < split.address.size(); ++i) {
                            int t = split.address[i];
                            Word evaluations = 0;
                            for (Word points = W; points; points &= points - 1) {
                                const int w = low(points);
                                if (parity(t & w)) evaluations |= Word(1) << w;
                            }
                            classes[evaluations] |= Word(1) << i;
                        }
                        assert(classes.size() == static_cast<unsigned>(std::popcount(W)));
                        auto& choices = split.activities[W];
                        for (const auto& [_, addresses] : classes) choices.push_back(addresses);
                    }
                }
                splits.push_back(std::move(split));
            }
        }
    }
};

struct Fourier {
    Word support = 0;
    std::array<int, 64> coefficients{};
    std::array<int, 64> signs{};
};

Fourier walsh(Word table, int size) {
    Fourier result;
    for (int t = 0; t < size; ++t) {
        result.signs[t] = 1 - 2 * ((table >> t) & 1);
        result.coefficients[t] = result.signs[t];
    }
    for (int step = 1; step < size; step <<= 1) {
        for (int start = 0; start < size; start += step << 1) {
            for (int p = start; p < start + step; ++p) {
                int a = result.coefficients[p], b = result.coefficients[p + step];
                result.coefficients[p] = a + b;
                result.coefficients[p + step] = a - b;
            }
        }
    }
    for (int t = 0; t < size; ++t) {
        if (result.coefficients[t]) result.support |= Word(1) << t;
    }
    return result;
}

struct Choices {
    std::vector<Word> masks;
    int size;
};

bool fit(std::vector<Choices> const& choices, unsigned i, int cap, std::array<int, 64>& loads) {
    if (i == choices.size()) return true;
    for (Word mask : choices[i].masks) {
        bool available = true;
        for (Word todo = mask; todo; todo &= todo - 1) {
            if (loads[low(todo)] >= cap) {
                available = false;
                break;
            }
        }
        if (!available) continue;
        for (Word todo = mask; todo; todo &= todo - 1) ++loads[low(todo)];
        if (fit(choices, i + 1, cap, loads)) return true;
        for (Word todo = mask; todo; todo &= todo - 1) --loads[low(todo)];
    }
    return false;
}

struct Value {
    int M;
    int L;
    int packing;
    int mass_lower;
};

Value evaluate(const Fourier& f, const Split& split, Count& checked, bool direct) {
    int M = 0;
    int mass = 0;
    std::array<int, 64> loads{};
    std::vector<Choices> choices;
    for (unsigned b = 0; b < split.bands.size(); ++b) {
        Word support = f.support & split.frequency_masks[b];
        if (!support && !direct) continue;
        Word active = 0;
        for (unsigned i = 0; i < split.address.size(); ++i) {
            int t = split.address[i];
            int from_spectrum = 0;
            for (Word todo = support; todo; todo &= todo - 1) {
                int frequency = low(todo);
                from_spectrum += parity(t & frequency) ? -f.coefficients[frequency] : f.coefficients[frequency];
            }
            if (direct) {
                int from_truth = 0;
                for (int u : split.physical_kernel)
                    from_truth += parity(split.bands[b] & u) ? -f.signs[t ^ u] : f.signs[t ^ u];
                assert(from_truth * static_cast<int>(split.address.size()) == from_spectrum);
                ++checked;
            }
            if (from_spectrum) {
                active |= Word(1) << i;
                ++loads[i];
            }
        }
        if (!support) continue;
        assert(active);
        Word fiber = span(translate(support, low(support)));
        assert((fiber & split.space) == fiber);
        M += std::popcount(fiber);
        Choices options{};
        for (Word piece : split.activities.at(fiber)) {
            if (piece & active) {
                assert((piece & active) == piece);
                options.masks.push_back(piece);
            }
        }
        assert(!options.masks.empty());
        options.size = std::popcount(options.masks.front());
        mass += options.size;
        choices.push_back(std::move(options));
    }
    int L = *std::max_element(loads.begin(), loads.end());
    const int address_size = split.address.size();
    int lower = (mass + address_size - 1) / address_size;
    assert(lower > 0 && lower <= L);
    std::sort(choices.begin(), choices.end(), [](const Choices& a, const Choices& b) {
        if (a.masks.size() != b.masks.size()) return a.masks.size() < b.masks.size();
        return a.size > b.size;
    });
    int optimum = lower;
    for (; optimum < L; ++optimum) {
        std::array<int, 64> selected{};
        if (fit(choices, 0, optimum, selected)) break;
    }
    return {M, L, optimum, lower};
}

struct Function {
    int support;
    int affine_rank;
    int best_raw;
    int best_pack;
    int balanced_raw;
    int balanced_pack;
    unsigned best_raw_count;
    unsigned best_pack_count;
    Word support_mask;
};

struct PhaseRecord {
    Word truth_min = 0;
    Word truth_max = 0;
    int minimum = 100;
    int maximum = 0;

    void check(Word truth, int value) {
        if (value < minimum || (value == minimum && truth < truth_min)) {
            minimum = value;
            truth_min = truth;
        }
        if (value > maximum || (value == maximum && truth < truth_max)) {
            maximum = value;
            truth_max = truth;
        }
    }

    void include(const PhaseRecord& other) {
        check(other.truth_min, other.minimum);
        check(other.truth_max, other.maximum);
    }
};

struct Run {
    Count tables = 0;
    Count tested = 0;
    Count direct = 0;
    Count pack_improves = 0;
    Count density_strict = 0;
    Count lower_insufficient = 0;
    std::map<std::array<int, 5>, Count> by_dim;
    std::map<std::array<int, 8>, Count> by_function;
    std::map<std::array<Word, 2>, std::array<PhaseRecord, 2>> phase;
    std::map<std::array<Word, 2>, std::array<PhaseRecord, 2>> exact_magnitude;
    std::map<std::array<int, 2>, Word> worst_balanced_fns;
};

Function evaluate_all(const Geometry& geom, Word truth, Run& result, bool direct, bool record_phases) {
    Fourier transform = walsh(truth, geom.size);
    const int s = std::popcount(transform.support);
    const int a = dimension(span(translate(transform.support, low(transform.support))));
    Function fn{s, a, 10000, 10000, 10000, 10000, 0, 0, transform.support};
    Word magnitudes = 0;
    if (record_phases) {
        assert(geom.size <= 16);
        for (int i = 0; i < geom.size; ++i)
            magnitudes |= Word(std::abs(transform.coefficients[i]) / 2) << (4 * i);
    }
    for (const Split& split : geom.splits) {
        Value value = evaluate(transform, split, result.direct, direct);
        const int raw = value.M * value.L;
        const int packed = value.M * value.packing;
        assert(1 <= value.packing && value.packing <= value.L && value.M >= s);
        assert(packed >= s);
        if (raw < fn.best_raw) fn.best_raw = raw, fn.best_raw_count = 0;
        if (packed < fn.best_pack) fn.best_pack = packed, fn.best_pack_count = 0;
        fn.best_raw_count += raw == fn.best_raw;
        fn.best_pack_count += packed == fn.best_pack;
        if (split.dim == geom.n / 2) {
            fn.balanced_raw = std::min(fn.balanced_raw, raw);
            fn.balanced_pack = std::min(fn.balanced_pack, packed);
        }
        ++result.tested;
        result.pack_improves += value.packing < value.L;
        result.density_strict += value.M > s;
        result.lower_insufficient += value.packing > value.mass_lower;
        ++result.by_dim[{split.dim, value.M, value.L, value.packing, value.mass_lower}];
        if (record_phases) {
            auto& phase = result.phase[{transform.support, split.space}];
            phase[0].check(truth, value.L);
            phase[1].check(truth, value.packing);
            auto& equal = result.exact_magnitude[{magnitudes, split.space}];
            equal[0].check(truth, value.L);
            equal[1].check(truth, value.packing);
        }
    }
    ++result.tables;
    ++result.by_function[{s, a, fn.best_raw, fn.best_pack, fn.balanced_raw, fn.balanced_pack,
                          int(fn.best_raw_count), int(fn.best_pack_count)}];
    auto key = std::array<int, 2>{s, fn.balanced_pack};
    auto [iterator, inserted] = result.worst_balanced_fns.try_emplace(key, truth);
    if (!inserted) iterator->second = std::min(iterator->second, truth);
    return fn;
}

template<unsigned N>
void profile_json(std::ostream& out, const std::map<std::array<int, N>, Count>& p) {
    out << "{";
    bool separator = false;
    for (const auto& [key, value] : p) {
        if (separator) out << ",";
        separator = true;
        out << "\"";
        for (unsigned j = 0; j < N; ++j) {
            if (j) out << ",";
            out << key[j];
        }
        out << "\":" << value;
    }
    out << "}";
}

void phases_json(std::ostream& out, std::map<std::array<Word, 2>, std::array<PhaseRecord, 2>> const& records) {
    out << "[";
    for (int kind = 0; kind < 2; ++kind) {
        if (kind) out << ",";
        int width = 0;
        Count varied = 0;
        std::array<Word, 2> witness{};
        for (const auto& [key, pair] : records) {
            int spread = pair[kind].maximum - pair[kind].minimum;
            varied += spread > 0;
            if (spread > width) {
                width = spread;
                witness = key;
            }
        }
        out << "{\"varying_same_support_and_kernel_cases\":" << varied << ",\"maximum_spread\":" << width;
        if (width) {
            const auto& record = records.at(witness)[kind];
            out << ",\"witness\":{\"first_key\":" << witness[0] << ",\"kernel\":" << witness[1]
                << ",\"value_min\":" << record.minimum << ",\"value_max\":" << record.maximum
                << ",\"truth_min\":" << record.truth_min << ",\"truth_max\":" << record.truth_max << "}";
        }
        out << "}";
    }
    out << "]";
}

void run_json(std::ostream& out, Run const& run) {
    out << "{\"truth_tables\":" << run.tables << ",\"kernel_tests\":" << run.tested
        << ",\"direct_physical_coefficients\":" << run.direct
        << ",\"strict_packing_improvements\":" << run.pack_improves
        << ",\"strict_affine_fill\":" << run.density_strict
        << ",\"strict_packing_over_volume_bound\":" << run.lower_insufficient
        << ",\"profile_dim_M_L_packing_mass_lower\":";
    profile_json<5>(out, run.by_dim);
    out << ",\"profile_s_affine_rank_best_raw_best_pack_balanced_raw_balanced_pack_winners_raw_winners_pack\":";
    profile_json<8>(out, run.by_function);
    out << ",\"phase_comparisons_raw_and_packed\":";
    phases_json(out, run.phase);
    out << ",\"exact_magnitude_comparisons_raw_and_packed\":";
    phases_json(out, run.exact_magnitude);
    out << ",\"smallest_truth_by_s_and_best_balanced_packed_product\":{";
    bool separator = false;
    for (const auto& [key, value] : run.worst_balanced_fns) {
        if (separator) out << ",";
        separator = true;
        out << "\"" << key[0] << "," << key[1] << "\":" << value;
    }
    out << "}}";
}

void merge(Run& into, Run const& other) {
    into.tables += other.tables;
    into.tested += other.tested;
    into.direct += other.direct;
    into.pack_improves += other.pack_improves;
    into.density_strict += other.density_strict;
    into.lower_insufficient += other.lower_insufficient;
    for (auto [key, count] : other.by_dim) into.by_dim[key] += count;
    for (auto [key, count] : other.by_function) into.by_function[key] += count;
    for (auto const& [key, pair] : other.phase) {
        auto& merged = into.phase[key];
        merged[0].include(pair[0]);
        merged[1].include(pair[1]);
    }
    for (const auto& [key, pair] : other.exact_magnitude) {
        auto& merged = into.exact_magnitude[key];
        merged[0].include(pair[0]);
        merged[1].include(pair[1]);
    }
    for (auto [key, truth] : other.worst_balanced_fns) {
        auto [iterator, inserted] = into.worst_balanced_fns.try_emplace(key, truth);
        if (!inserted) iterator->second = std::min(iterator->second, truth);
    }
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) return 2;
    const int workers = std::stoi(argv[1]);
    const int samples = std::stoi(argv[3]);
    if (workers <= 0 || samples < 0) return 2;
    std::ofstream out(argv[2]);
    out << "{\"schema\":\"a28-r072-all-kernel-minimization-v1\",\"workers\":" << workers
        << ",\"exhaustive\":[\n";
    for (int n = 1; n <= 4; ++n) {
        auto start = std::chrono::steady_clock::now();
        Geometry g(n);
        const Word tables = Word(1) << g.size;
        std::vector<Run> independent(workers);
        std::vector<std::thread> active;
        for (int j = 0; j < workers; ++j) {
            active.emplace_back([&, j] {
                for (Word truth = tables * j / workers; truth < tables * (j + 1) / workers; ++truth) {
                    Function value = evaluate_all(g, truth, independent[j], true, true);
                    assert(value.best_raw == value.support && value.best_pack == value.support);
                }
            });
        }
        for (auto& process : active) process.join();
        Run summary;
        for (const auto& run : independent) merge(summary, run);
        if (n > 1) out << ",\n";
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        out << "{\"n\":" << n << ",\"kernels\":" << g.splits.size()
            << ",\"seconds\":" << elapsed << ",\"summary\":";
        run_json(out, summary);
        out << "}";
        std::cout << "exhaustive " << n << ": " << summary.tested << " " << elapsed << std::endl;
    }
    out << "]";
    if (samples || argc == 5) {
        const int n = 5;
        auto start = std::chrono::steady_clock::now();
        Geometry g(n);
        constexpr Word seed = 720718;
        std::mt19937_64 rng(seed);
        std::vector<Word> tables;
        for (int j = 0; j < samples; ++j) tables.push_back(rng() & ((Word(1) << g.size) - 1));
        if (argc == 5) {
            std::ifstream extra(argv[4]);
            Word mask;
            while (extra >> mask) tables.push_back(mask);
            assert(extra.eof());
        }
        const int total = tables.size();
        std::vector<Run> independent(workers);
        std::vector<std::thread> active;
        std::vector<Function> values(total);
        for (int j = 0; j < workers; ++j) {
            active.emplace_back([&, j] {
                for (int i = total * j / workers; i < total * (j + 1) / workers; ++i)
                    values[i] = evaluate_all(g, tables[i], independent[j], false, false);
            });
        }
        for (auto& process : active) process.join();
        Run summary;
        for (const auto& run : independent) merge(summary, run);
        unsigned raw_index = 0, pack_index = 0, earliest_index = 0;
        Count obstructions = 0;
        Count additional_obstructions = 0;
        for (unsigned i = 0; i < values.size(); ++i) {
            const auto& value = values[i];
            if (value.best_pack != value.support) {
                ++obstructions;
                if (i >= static_cast<unsigned>(samples)) ++additional_obstructions;
                if (values[earliest_index].best_pack == values[earliest_index].support || tables[i] < tables[earliest_index])
                    earliest_index = i;
            }
            if (value.best_raw * values[raw_index].support > value.support * values[raw_index].best_raw
                || (value.best_raw * values[raw_index].support == value.support * values[raw_index].best_raw
                    && tables[i] < tables[raw_index])) raw_index = i;
            if (value.best_pack * values[pack_index].support > value.support * values[pack_index].best_pack
                || (value.best_pack * values[pack_index].support == value.support * values[pack_index].best_pack
                    && tables[i] < tables[pack_index])) pack_index = i;
        }
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        out << ",\"n5_search\":{\"random_samples\":" << samples << ",\"additional_tables\":" << (total - samples)
            << ",\"seed\":" << seed
            << ",\"packed_obstructions\":" << obstructions
            << ",\"additional_packed_obstructions\":" << additional_obstructions
            << ",\"seconds\":" << elapsed << ",\"summary\":";
        run_json(out, summary);
        out << ",\"maximum_raw_maximum_packed_and_smallest_obstruction_witnesses\":[";
        bool separator = false;
        for (unsigned i : {raw_index, pack_index, earliest_index}) {
            if (separator) out << ",";
            separator = true;
            auto fn = values[i];
            out << "{\"truth\":" << tables[i] << ",\"s\":" << fn.support << ",\"support_mask\":" << fn.support_mask
                << ",\"best_raw\":" << fn.best_raw << ",\"best_packed\":" << fn.best_pack
                << ",\"balanced_raw\":" << fn.balanced_raw << ",\"balanced_packed\":" << fn.balanced_pack << "}";
        }
        out << "]}";
        std::cout << "n5 " << total << " " << elapsed << std::endl;
    }
    out << "}\n";
}
