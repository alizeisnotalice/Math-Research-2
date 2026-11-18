#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <thread>
#include <vector>

using Mask = unsigned;
using Count = unsigned long long;

int parity(int x) { return std::popcount(static_cast<unsigned>(x)) & 1; }
int first(Mask mask) { return std::countr_zero(mask); }

Mask shift(Mask set, int offset) {
    Mask out = 0;
    for (; set; set &= set - 1) out |= 1u << (first(set) ^ offset);
    return out;
}

Mask adjoin(Mask space, int vector) { return space | shift(space, vector); }

int rank(Mask points) {
    Mask space = 1;
    for (; points; points &= points - 1) space = adjoin(space, first(points));
    return std::countr_zero(static_cast<unsigned>(std::popcount(space)));
}

struct Plane {
    Mask space;
    int dimension;
    std::array<int, 16> reduction{};
    std::vector<int> normals;
    std::vector<int> shifts;
    std::vector<int> contained;
    std::vector<int> representatives;
};

struct Core {
    Mask affine;
    int offset;
    int direction;
};

struct Geometry {
    int n;
    int size;
    std::vector<Plane> planes;
    std::vector<Core> cores;

    explicit Geometry(int n) : n(n), size(1 << n) {
        for (Mask s = 1; s < (1u << size); s += 2) {
            const int count = std::popcount(s);
            if (!std::has_single_bit(static_cast<unsigned>(count)) || rank(s) != std::countr_zero(static_cast<unsigned>(count))) continue;
            Plane p{};
            p.space = s;
            p.dimension = std::countr_zero(static_cast<unsigned>(count));
            for (int i = 0; i < size; ++i) {
                p.reduction[i] = first(shift(s, i));
                bool orthogonal = true;
                for (Mask left = s; left; left &= left - 1) orthogonal &= !parity(i & first(left));
                if (orthogonal) p.normals.push_back(i);
                if (p.reduction[i] == i) p.representatives.push_back(i);
            }
            Mask coset = 0;
            for (int i : p.normals) coset |= 1u << i;
            for (int i = 0; i < size; ++i) if (first(shift(coset, i)) == i) p.shifts.push_back(i);
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) assert(p.reduction[i ^ j] == (p.reduction[i] ^ p.reduction[j]));
            }
            planes.push_back(std::move(p));
        }
        for (int i = 0; i < static_cast<int>(planes.size()); ++i) {
            for (int j = 0; j < static_cast<int>(planes.size()); ++j) {
                if ((planes[i].space & planes[j].space) == planes[j].space) planes[i].contained.push_back(j);
            }
            for (int offset : planes[i].representatives) cores.push_back({shift(planes[i].space, offset), offset, i});
        }
    }
};

struct Combined {
    Mask support;
    Mask odd;
    std::vector<std::array<Mask, 16>> restrictions;
};

Combined exact_restrictions(const Geometry& g, unsigned truth, Count& checks) {
    std::array<int, 16> signs{}, spectrum{};
    for (int x = 0; x < g.size; ++x) signs[x] = 1 - 2 * ((truth >> x) & 1);
    int divisor = 0;
    Mask support = 0;
    for (int frequency = 0; frequency < g.size; ++frequency) {
        for (int x = 0; x < g.size; ++x) spectrum[frequency] += parity(x & frequency) ? -signs[x] : signs[x];
        divisor = std::gcd(divisor, spectrum[frequency]);
        if (spectrum[frequency]) support |= 1u << frequency;
    }
    Mask odd = 0;
    for (int frequency = 0; frequency < g.size; ++frequency) {
        if (spectrum[frequency] / divisor & 1) odd |= 1u << frequency;
    }
    Combined full{support, odd, {}};
    for (const Plane& p : g.planes) {
        std::array<Mask, 16> restrictions{};
        for (int t : p.shifts) {
            for (int freq : p.representatives) {
                int physical = 0;
                for (int h : p.normals) physical += parity(freq & h) ? -signs[t ^ h] : signs[t ^ h];
                int folded = 0;
                for (Mask left = p.space; left; left &= left - 1) {
                    int point = freq ^ first(left);
                    folded += parity(t & point) ? -spectrum[point] : spectrum[point];
                }
                assert(physical * std::popcount(p.space) == folded);
                if (physical) restrictions[t] |= 1u << freq;
                ++checks;
            }
        }
        full.restrictions.push_back(restrictions);
    }
    return full;
}

struct Witness {
    unsigned truth = 0;
    Mask core = 0;
    Mask differences = 0;
    int r = 0;
    int q = 0;
    int cost = 0;
    int best_rank = 0;
    int best_support = 0;
};

struct Profile {
    Count cores = 0;
    Count nonempty_tails = 0;
    Count kernels = 0;
    Count affine_restrictions = 0;
    Count collision_free_kernels = 0;
    Count lower_loss_at_h0 = 0;
    Count strict_per_coset_loss = 0;
    Count no_full_rank_at_h0 = 0;
    Count deduplicated_budget_failures = 0;
    Count small_compression_nontrivial = 0;
    Count zero_loss_nontrivial = 0;
    std::map<std::array<int, 4>, Count> shapes;
    Witness example;
    Witness deduplicated_example;
};

void earliest(Witness& into, const Witness& from) {
    if (from.core && (!into.core || from.truth < into.truth ||
                      (from.truth == into.truth && from.core < into.core)))
        into = from;
}

void merge(Profile& into, const Profile& from) {
    into.cores += from.cores;
    into.nonempty_tails += from.nonempty_tails;
    into.kernels += from.kernels;
    into.affine_restrictions += from.affine_restrictions;
    into.collision_free_kernels += from.collision_free_kernels;
    into.lower_loss_at_h0 += from.lower_loss_at_h0;
    into.strict_per_coset_loss += from.strict_per_coset_loss;
    into.no_full_rank_at_h0 += from.no_full_rank_at_h0;
    into.deduplicated_budget_failures += from.deduplicated_budget_failures;
    into.small_compression_nontrivial += from.small_compression_nontrivial;
    into.zero_loss_nontrivial += from.zero_loss_nontrivial;
    for (auto [shape, count] : from.shapes) into.shapes[shape] += count;
    earliest(into.example, from.example);
    earliest(into.deduplicated_example, from.deduplicated_example);
}

void evaluate(const Geometry& geom, unsigned truth, const Combined& full, const Core& core, Profile& stats) {
    const Plane& parent = geom.planes[core.direction];
    const Mask shifted = shift(full.support, core.offset);
    const Mask tail = shifted & ~parent.space;
    const int q = std::popcount(tail);
    const int core_dimension = parent.dimension;
    const int M = 1 << core_dimension;
    std::array<Mask, 16> bands{};
    for (Mask current = tail; current; current &= current - 1) {
        const int point = first(current);
        bands[parent.reduction[point]] |= 1u << point;
    }
    Mask external = 0;
    for (int gamma : parent.representatives) if (bands[gamma]) external |= 1u << gamma;
    const int r = rank(external);
    Mask basis_span = 1;
    Mask basis_labels = 0;
    Mask differences = 0;
    int B = 0;
    for (Mask current = external; current; current &= current - 1) {
        const int gamma = first(current);
        if ((basis_span >> gamma) & 1u) continue;
        basis_span = adjoin(basis_span, gamma);
        basis_labels |= 1u << gamma;
        const int representative = first(bands[gamma]);
        differences |= shift(bands[gamma], representative) & ~1u;
        B += std::popcount(bands[gamma]) - 1;
    }
    const int D = std::popcount(differences);
    assert((differences & ~parent.space) == 0);
    const int h_target = std::bit_width(static_cast<unsigned>(D + 1)) - 1;
    ++stats.cores;
    if (q) {
        ++stats.nonempty_tails;
        ++stats.shapes[{M, q, r, D}];
    }
    std::array<Count, 5> kernels{}, affine_counts{}, marked_losses{}, genuine_losses{};
    bool clean = false;
    int best_at_h0 = 0;
    int max_support_at_h0 = 0;
    for (int index : parent.contained) {
        const Plane& p = geom.planes[index];
        const int h = core_dimension - p.dimension;
        const int L = 1 << h;
        ++stats.kernels;
        ++kernels[h];
        if (!(p.space & differences)) {
            ++stats.collision_free_kernels;
            if (h == h_target) clean = true;
        }
        int marked_collisions = 0;
        for (Mask current = basis_labels; current; current &= current - 1) {
            const int gamma = first(current);
            marked_collisions += ((shift(bands[gamma], first(bands[gamma])) & ~1u) & p.space) != 0;
        }
        marked_losses[h] += marked_collisions;
        for (int t : p.shifts) {
            const Mask restricted = shift(full.restrictions[index][t], p.reduction[core.offset]);
            Mask image = 0;
            for (Mask current = restricted; current; current &= current - 1)
                image |= 1u << parent.reduction[first(current)];
            image &= ~1u;
            assert((image & ~external) == 0);
            const int rank_kept = rank(image);
            const int loss = r - rank_kept;
            assert(loss >= 0 && loss <= marked_collisions);
            assert(std::popcount(restricted) <= q + L);
            if (!(p.space & differences)) assert(loss == 0);
            if (loss) ++stats.strict_per_coset_loss;
            genuine_losses[h] += loss;
            ++affine_counts[h];
            ++stats.affine_restrictions;
            if (!h) {
                if (rank_kept > best_at_h0) {
                    best_at_h0 = rank_kept;
                    max_support_at_h0 = std::popcount(restricted);
                } else if (rank_kept == best_at_h0) {
                    max_support_at_h0 = std::max(max_support_at_h0, std::popcount(restricted));
                }
            }
        }
    }
    assert(clean);
    if (q && h_target < core_dimension) ++stats.zero_loss_nontrivial;
    for (int h = 0; h <= core_dimension; ++h) {
        const int L = 1 << h;
        const int K = M / L;
        assert(genuine_losses[h] * L <= affine_counts[h] * B);
        assert(marked_losses[h] * (M - 1) <= kernels[h] * B * (K - 1));
        if (M > 1) {
            const Count count = affine_counts[h];
            const Count loss = genuine_losses[h];
            assert(2 * loss * (M - 1) <= count * B * (K - 1));
            assert((count * r - loss) * (r * (M - 1) + B * (K - 1))
                   >= count * r * r * (M - 1));
        }
        if (genuine_losses[h] * L > affine_counts[h] * D) {
            ++stats.deduplicated_budget_failures;
            if (h == 0)
                earliest(stats.deduplicated_example,
                         {truth, core.affine, differences, r, q, B, best_at_h0, max_support_at_h0});
        }
        if (h < h_target && B < r * L && q) ++stats.small_compression_nontrivial;
    }
    if (genuine_losses[0]) ++stats.lower_loss_at_h0;
    if (best_at_h0 < r) {
        ++stats.no_full_rank_at_h0;
        earliest(stats.example, {truth, core.affine, differences, r, q, B, best_at_h0, max_support_at_h0});
    }
}

struct Run {
    Count tables = 0;
    Count checks = 0;
    Profile all;
    Profile highest;
};

Run run_range(const Geometry& g, unsigned begin, unsigned end) {
    Run result;
    for (unsigned truth = begin; truth < end; ++truth) {
        ++result.tables;
        Combined full = exact_restrictions(g, truth, result.checks);
        bool odd_found = false;
        for (const Core& c : g.cores) {
            if ((c.affine & full.support) != c.affine) continue;
            evaluate(g, truth, full, c, result.all);
            if (c.affine == full.odd) {
                evaluate(g, truth, full, c, result.highest);
                odd_found = true;
            }
        }
        assert(odd_found);
    }
    return result;
}

void witness_json(std::ostream& out, const Witness& example, const char* name) {
    if (!example.core) return;
    out << ",\"" << name << "\":{\"truth_mask\":" << example.truth
        << ",\"core_point_mask\":" << example.core << ",\"difference_point_mask\":" << example.differences
        << ",\"tail_rank\":" << example.r << ",\"tail_points\":" << example.q
        << ",\"collision_budget\":" << example.cost << ",\"best_external_rank_at_h0\":" << example.best_rank
        << ",\"max_support_among_best_at_h0\":" << example.best_support << "}";
}

void profile_json(std::ostream& out, const Profile& p) {
    out << "{\"cores\":" << p.cores << ",\"nonempty_tails\":" << p.nonempty_tails
        << ",\"kernels\":" << p.kernels << ",\"affine_restrictions\":" << p.affine_restrictions
        << ",\"collision_free_kernels\":" << p.collision_free_kernels
        << ",\"cores_with_loss_at_h0\":" << p.lower_loss_at_h0
        << ",\"individual_cosets_with_loss\":" << p.strict_per_coset_loss
        << ",\"cores_with_no_full_rank_at_h0\":" << p.no_full_rank_at_h0
        << ",\"deduplicated_budget_failures\":" << p.deduplicated_budget_failures
        << ",\"nontrivial_positive_expected_rank_at_h_below_target\":" << p.small_compression_nontrivial
        << ",\"can_reduce_core_preserving_full_external_rank\":" << p.zero_loss_nontrivial;
    witness_json(out, p.example, "first_full_collapse_obstruction");
    witness_json(out, p.deduplicated_example, "first_deduplicated_budget_failure_at_h0");
    out << ",\"shapes_M_q_r_D\":{";
    bool separator = false;
    for (auto [shape, count] : p.shapes) {
        if (separator) out << ",";
        separator = true;
        out << "\"" << shape[0] << "," << shape[1] << "," << shape[2] << "," << shape[3] << "\":" << count;
    }
    out << "}}";
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int workers = std::stoi(argv[1]);
    if (workers < 1) return 2;
    std::ofstream out(argv[2]);
    out << "{\"schema\":\"a28-r063-exact-kernel-compression-v1\",\"workers\":" << workers << ",\"dimensions\":[\n";
    for (int n = 1; n <= 4; ++n) {
        auto begin = std::chrono::steady_clock::now();
        const Geometry g(n);
        const unsigned total = 1u << g.size;
        std::vector<Run> partial(workers);
        std::vector<std::thread> threads;
        for (int i = 0; i < workers; ++i)
            threads.emplace_back([&, i] { partial[i] = run_range(g, total * i / workers, total * (i + 1) / workers); });
        for (auto& thread : threads) thread.join();
        Run combined;
        for (const Run& run : partial) {
            combined.tables += run.tables;
            combined.checks += run.checks;
            merge(combined.all, run.all);
            merge(combined.highest, run.highest);
        }
        assert(combined.tables == total);
        if (n > 1) out << ",\n";
        out << "{\"n\":" << n << ",\"truth_tables\":" << total << ",\"affine_flat_catalog_size\":" << g.cores.size()
            << ",\"direct_restriction_coefficients\":" << combined.checks << ",\"all_full_affine_subcores\":";
        profile_json(out, combined.all);
        out << ",\"highest_odd_layer_cores\":";
        profile_json(out, combined.highest);
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        out << ",\"elapsed_seconds\":" << elapsed << "}";
        out.flush();
        std::cout << n << " " << total << " " << combined.all.cores << " " << elapsed << std::endl;
    }
    out << "\n]}\n";
}
