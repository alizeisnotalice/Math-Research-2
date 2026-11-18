#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using Word = std::uint64_t;

constexpr int starts[4] = {0, 2, 3, 5};
constexpr int ends[4] = {1, 4, 7, 6};
constexpr const char* labels[5] = {
    "address", "local_ancestor", "within_line_siblings",
    "within_line_nonsiblings", "across_lines"
};

void require(bool condition, const char* explanation) {
    if (!condition) {
        throw std::runtime_error(explanation);
    }
}

int parity(Word value) {
    return std::popcount(value) & 1;
}

struct LinearRank {
    std::array<Word, 64> rows{};
    int value = 0;

    void insert(Word point) {
        while (point) {
            const int pivot = std::bit_width(point) - 1;
            if (rows[pivot]) {
                point ^= rows[pivot];
            } else {
                rows[pivot] = point;
                ++value;
                break;
            }
        }
    }
};

using Spectrum = std::vector<std::pair<Word, int>>;

Spectrum ancestor_terms(int node, int offset) {
    Spectrum terms{{0, 1}};
    while (node) {
        const int parent = (node - 1) / 2;
        const int branch_sign = node == 2 * parent + 1 ? 1 : -1;
        const std::size_t old_size = terms.size();
        for (std::size_t i = 0; i < old_size; ++i) {
            terms.emplace_back(terms[i].first | (Word{1} << (offset + parent)),
                               terms[i].second * branch_sign);
        }
        node = parent;
    }
    return terms;
}

Spectrum build_parent(int height) {
    const int block_nodes = (1 << height) - 1;
    const int bottom = (1 << (height - 1)) - 1;
    Spectrum spectrum;
    for (int group = 0; group < 4; ++group) {
        const int offset = 3 + group * block_nodes;
        for (int node = bottom; node < block_nodes; ++node) {
            const Word terminal = Word{1} << (offset + node);
            for (const auto& [ancestors, sign] : ancestor_terms(node, offset)) {
                for (Word selector = 0; selector < 8; ++selector) {
                    if (parity(selector & (starts[group] ^ ends[group])) == 0) {
                        const int selector_sign = 1 - 2 * parity(selector & starts[group]);
                        spectrum.emplace_back(terminal | ancestors | selector,
                                              sign * selector_sign);
                    }
                }
            }
        }
    }
    return spectrum;
}

int direct_value(Word point, int height) {
    const int selector = point & 7;
    int group = 0;
    while (selector != starts[group] && selector != ends[group]) {
        ++group;
        require(group < 4, "physical selector partition");
    }
    const int block_nodes = (1 << height) - 1;
    const int offset = 3 + group * block_nodes;
    const int bottom = (1 << (height - 1)) - 1;
    int node = 0;
    while (node < bottom) {
        node = 2 * node + 1 + ((point >> (offset + node)) & 1);
    }
    return 1 - 2 * ((point >> (offset + node)) & 1);
}

void check_direct(const std::unordered_map<Word, int>& parent, int height, int dimension, int norm) {
    std::vector<int> values(1 << dimension);
    for (std::size_t point = 0; point < values.size(); ++point) {
        values[point] = direct_value(point, height);
    }
    for (std::size_t width = 1; width < values.size(); width *= 2) {
        for (std::size_t start = 0; start < values.size(); start += 2 * width) {
            for (std::size_t i = start; i < start + width; ++i) {
                const int left = values[i], right = values[i + width];
                values[i] = left + right;
                values[i + width] = left - right;
            }
        }
    }
    for (std::size_t frequency = 0; frequency < values.size(); ++frequency) {
        const auto found = parent.find(frequency);
        const int coefficient = found == parent.end() ? 0 : found->second;
        require(values[frequency] * norm == coefficient * static_cast<int>(values.size()),
                "complete direct parent Walsh");
    }
}

struct Profile {
    int kind;
    int doubles;
    int child_sparsity;
    int intersection;
    int squared_slack;

    auto operator<=>(const Profile&) const = default;
};

Profile evaluate_cut(const Spectrum& spectrum,
                     const std::unordered_map<Word, int>& parent,
                     Word query, Word terminal_mask,
                     int height, int dimension, int norm) {
    const Word pivot = Word{1} << (std::bit_width(query) - 1);
    const Word flag = Word{1} << dimension;
    LinearRank zero, one, union_span;
    int size_zero = 0, size_one = 0, doubles = 0;
    for (const auto& [point, sign] : spectrum) {
        const auto found = parent.find(point ^ query);
        const bool paired = found != parent.end();
        if (paired && (point & pivot)) {
            continue;
        }
        const Word representative = point ^ ((point & pivot) ? query : 0);
        const Word augmented = representative | flag;
        union_span.insert(augmented);
        const int other_sign = paired ? found->second : 0;
        if (sign + other_sign) {
            zero.insert(augmented);
            ++size_zero;
        }
        if (sign - other_sign) {
            one.insert(augmented);
            ++size_one;
        }
        doubles += paired;
    }
    require(doubles > 0 && doubles % 2 == 0, "strict cut and balanced signs");
    require(size_zero == size_one && size_one == norm * norm - 3 * doubles / 2,
            "flat cut exact child sparsities");
    require(union_span.value == 2 * norm - 2, "rank of projected parent");
    const int intersection = zero.value + one.value - union_span.value;
    const int squared_slack = (intersection + 2 * norm) * (intersection + 2 * norm)
                            - 16 * size_one;
    require(squared_slack < 0, "no C=2 payable cut");

    const int block_nodes = (1 << height) - 1;
    const Word terminals = query & terminal_mask;
    if (!terminals && query < 8) {
        require(intersection == norm && doubles == norm * norm / 4,
                "pure selector formula");
        return {0, doubles, size_one, intersection, squared_slack};
    }
    if (!terminals) {
        const int largest_bit = std::bit_width(query) - 1;
        const int local = (largest_bit - 3) % block_nodes;
        const int depth = std::bit_width(static_cast<unsigned>(local + 1)) - 1;
        const int affected = 1 << (height - 1 - depth);
        require(intersection == 2 * norm - 2 * affected && doubles == norm * affected / 2,
                "local ancestor formula");
        return {1, doubles, size_one, intersection, squared_slack};
    }
    require(std::popcount(terminals) == 2, "two terminal labels");
    const int first = std::countr_zero(terminals) - 3;
    const int second = std::bit_width(terminals) - 1 - 3;
    if (first / block_nodes != second / block_nodes) {
        require(intersection == 2 * norm - 2 && doubles == 2,
                "across lines formula");
        return {4, doubles, size_one, intersection, squared_slack};
    }
    if ((first % block_nodes - 1) / 2 == (second % block_nodes - 1) / 2) {
        require(intersection == 2 * norm - 4 && doubles == norm,
                "within-line sibling formula");
        return {2, doubles, size_one, intersection, squared_slack};
    }
    require(intersection == 2 * norm - 2 && doubles <= norm / 2,
            "within-line non-sibling formula");
    return {3, doubles, size_one, intersection, squared_slack};
}

void audit_height(int height) {
    const auto started = std::chrono::steady_clock::now();
    const int block_nodes = (1 << height) - 1;
    const int bottom = (1 << (height - 1)) - 1;
    const int dimension = 3 + 4 * block_nodes;
    const int norm = 4 * (1 << (height - 1));
    const Spectrum spectrum = build_parent(height);
    require(spectrum.size() == static_cast<std::size_t>(norm * norm), "parent sparsity");
    std::unordered_map<Word, int> parent;
    LinearRank lifted;
    for (const auto& [point, sign] : spectrum) {
        require(parent.emplace(point, sign).second, "parent bands disjoint");
        lifted.insert(point | (Word{1} << dimension));
    }
    require(lifted.value == 2 * norm - 1, "augmented parent rank");
    if (dimension <= 15) {
        check_direct(parent, height, dimension, norm);
    }
    std::unordered_map<Word, int> correlations;
    correlations.reserve(spectrum.size() * (spectrum.size() - 1) / 2);
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            correlations[spectrum[i].first ^ spectrum[j].first]
                += spectrum[i].second * spectrum[j].second;
        }
    }
    Word terminal_mask = 0;
    for (int group = 0; group < 4; ++group) {
        for (int node = bottom; node < block_nodes; ++node) {
            terminal_mask |= Word{1} << (3 + group * block_nodes + node);
        }
    }
    std::map<Profile, int> census;
    std::map<int, int> class_max_h;
    int maximum_slack = std::numeric_limits<int>::min();
    for (const auto& [query, correlation] : correlations) {
        require(query && correlation == 0, "all nonconstant coefficients of F squared vanish");
        const auto profile = evaluate_cut(spectrum, parent, query, terminal_mask,
                                          height, dimension, norm);
        ++census[profile];
        class_max_h[profile.kind] = std::max(class_max_h[profile.kind], profile.intersection);
        maximum_slack = std::max(maximum_slack, profile.squared_slack);
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\"height\":" << height << ",\"ambient_dimension\":" << dimension
              << ",\"N\":" << norm << ",\"sparsity\":" << spectrum.size()
              << ",\"affine_rank\":" << lifted.value - 1
              << ",\"strict_cuts\":" << correlations.size()
              << ",\"full_truth_table_walsh\":" << (dimension <= 15 ? "true" : "false")
              << ",\"exact_boolean_squared_spectrum\":true,\"maximum_C2_squared_slack\":"
              << maximum_slack << ",\"class_max_h\":{";
    bool first = true;
    for (const auto& [kind, maximum] : class_max_h) {
        std::cout << (first ? "" : ",") << "\"" << labels[kind] << "\":" << maximum;
        first = false;
    }
    std::cout << "},\"profiles\":[";
    first = true;
    for (const auto& [p, count] : census) {
        std::cout << (first ? "" : ",") << "{\"type\":\"" << labels[p.kind]
                  << "\",\"P\":" << p.doubles
                  << ",\"child_sparsity\":" << p.child_sparsity
                  << ",\"h\":" << p.intersection << ",\"squared_slack\":"
                  << p.squared_slack << ",\"count\":" << count << "}";
        first = false;
    }
    std::cout << "],\"elapsed_seconds\":" << elapsed << "}";
}

int main() {
    std::cout << "{\"schema\":\"a28-r049-paired-lines-cuts-v1\",\"compiler\":\"" << __VERSION__
              << "\",\"arithmetic\":\"exact integer, 64-bit GF(2)\",\"heights\":[";
    for (int height = 1; height <= 4; ++height) {
        if (height > 1) {
            std::cout << ",\n";
        }
        audit_height(height);
    }
    std::cout << "]}\n";
}
