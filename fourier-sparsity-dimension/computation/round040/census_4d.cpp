#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <map>
#include <thread>
#include <vector>

using Mask = uint16_t;
using Code = std::vector<Mask>;
using Hist = std::array<uint64_t, 17>;

Code span(const Code& generators) {
    Code out{0};
    for (Mask g : generators) {
        size_t before = out.size();
        for (size_t i = 0; i < before; ++i) out.push_back(out[i] ^ g);
    }
    std::sort(out.begin(), out.end());
    return out;
}

Code basis(const Code& space) {
    Code out, generated{0};
    for (Mask g : space) {
        if (!std::binary_search(generated.begin(), generated.end(), g)) {
            out.push_back(g);
            generated = span(out);
        }
    }
    return out;
}

std::vector<std::vector<Code>> subspaces(const Code& admissible, int depth) {
    std::vector<std::vector<Code>> result(depth + 1);
    result[0].push_back(Code{0});
    for (int d = 1; d <= depth; ++d) {
        std::set<Code> next;
        for (const auto& small : result[d - 1]) {
            for (Mask v : admissible) {
                if (std::binary_search(small.begin(), small.end(), v)) continue;
                Code candidate = small;
                bool legal = true;
                for (Mask x : small) {
                    Mask y = x ^ v;
                    if (!std::binary_search(admissible.begin(), admissible.end(), y)) {
                        legal = false;
                        break;
                    }
                    candidate.push_back(y);
                }
                if (legal) {
                    std::sort(candidate.begin(), candidate.end());
                    next.insert(candidate);
                }
            }
        }
        result[d].assign(next.begin(), next.end());
    }
    return result;
}

std::array<Mask, 16> characters() {
    std::array<Mask, 16> chars{};
    for (unsigned y = 0; y < 16; ++y) {
        for (unsigned x = 0; x < 16; ++x) {
            if (std::popcount(x & y) % 2) chars[y] |= Mask(1u << x);
        }
    }
    return chars;
}

bool bent(Mask phase, const std::array<Mask, 16>& chars) {
    for (Mask c : chars) {
        int w = std::popcount(Mask(phase ^ c));
        if (w != 6 && w != 10) return false;
    }
    return true;
}

std::array<Mask, 16> point_action(unsigned v, bool translation) {
    std::array<Mask, 16> move{};
    for (unsigned x = 0; x < 16; ++x) {
        unsigned y = translation
            ? x ^ v
            : x ^ (v * (std::popcount((((v & 5u) << 1) | ((v & 10u) >> 1)) & x) % 2));
        move[x] = Mask(1u << y);
    }
    return move;
}

Code transform(const Code& code, const std::array<Mask, 16>& move) {
    Code out;
    for (Mask g : code) {
        Mask h = 0;
        for (unsigned x = 0; x < 16; ++x) {
            if ((g >> x) & 1u) h |= move[x];
        }
        out.push_back(h);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::map<int, int> orbit_sizes(const std::set<Code>& codes) {
    std::vector<std::array<Mask, 16>> moves;
    for (unsigned v = 1; v < 16; ++v) moves.push_back(point_action(v, false));
    for (unsigned v : {1, 2, 4, 8}) moves.push_back(point_action(v, true));
    std::set<Code> unseen = codes;
    std::map<int, int> sizes;
    while (!unseen.empty()) {
        std::vector<Code> orbit{*unseen.begin()};
        unseen.erase(unseen.begin());
        for (size_t i = 0; i < orbit.size(); ++i) {
            for (const auto& move : moves) {
                auto image = transform(orbit[i], move);
                if (!codes.count(image)) std::abort();
                if (unseen.erase(image)) orbit.push_back(std::move(image));
            }
        }
        ++sizes[int(orbit.size())];
    }
    return sizes;
}

struct Kernel {
    Code basis;
    Code coset_representatives;
    Mask support = 0;
};

struct Result {
    std::array<Hist, 8> census{};
    std::array<Hist, 8> with_e4{};
    std::array<Code, 8> witness;
    unsigned legal = 0;
    std::array<unsigned, 3> projections{};
    std::map<int, int> d4_orbits;
    std::map<int, int> d5_orbits;
};

bool has_e4(const Code& cube, const Code& hyperplanes) {
    for (Mask h : hyperplanes) {
        if (!std::binary_search(cube.begin(), cube.end(), h)) continue;
        unsigned local = 0;
        for (Mask v : cube) if ((v & ~h) == 0) ++local;
        if (local == 8) return true;
    }
    return false;
}

Result run_seed(Mask f, const Code& quad, const Code& affine,
                const std::array<Mask, 16>& chars,
                const std::vector<std::vector<Kernel>>& kernels,
                const Code& hyperplanes) {
    Result result;
    std::array<bool, 65536> legal{};
    for (unsigned table = 0; table < 65536; ++table) {
        legal[table] = bent(Mask(table) ^ f, chars);
        result.legal += legal[table];
    }
    Code good;
    for (Mask q : quad) if (legal[q]) good.push_back(q);
    if (result.legal != good.size() * affine.size()) std::abort();
    for (Mask q : good) for (Mask a : affine) if (!legal[q ^ a]) std::abort();

    auto projected = subspaces(good, 3);
    if (!projected[3].empty()) std::abort();
    std::set<Code> d4_min, d5_min;
    for (int r = 0; r <= 2; ++r) {
        result.projections[r] = unsigned(projected[r].size());
        for (const auto& projection : projected[r]) {
            auto qb = basis(projection);
            for (int t = 0; t <= 5; ++t) {
                int d = t + r;
                if (d < 4) continue;
                for (const auto& k : kernels[t]) {
                    Code empty{0};
                    const auto& first = r >= 1 ? k.coset_representatives : empty;
                    const auto& second = r >= 2 ? k.coset_representatives : empty;
                    for (Mask a1 : first) {
                        Mask w1 = r >= 1 ? qb[0] ^ a1 : 0;
                        for (Mask a2 : second) {
                            Mask w2 = r >= 2 ? qb[1] ^ a2 : 0;
                            Mask support = k.support | w1 | w2;
                            unsigned n = std::popcount(support);
                            Code generators = k.basis;
                            if (r >= 1) generators.push_back(w1);
                            if (r >= 2) generators.push_back(w2);
                            auto cube = span(generators);
                            Mask actual_support = 0;
                            for (Mask w : cube) {
                                if (!legal[w]) std::abort();
                                actual_support |= w;
                            }
                            if (actual_support != support) std::abort();
                            bool e4 = (d == 4 || d == 5) && has_e4(cube, hyperplanes);
                            ++result.census[d][n];
                            result.with_e4[d][n] += e4;
                            if (result.witness[d].empty() || n < unsigned(result.witness[d][0])) {
                                result.witness[d] = Code{Mask(n)};
                                result.witness[d].insert(result.witness[d].end(), generators.begin(), generators.end());
                            }
                            if (f == 0x7888 && d == 4 && n == 12) d4_min.insert(cube);
                            if (f == 0x7888 && d == 5 && n <= 15) d5_min.insert(cube);
                        }
                    }
                }
            }
        }
    }
    if (f == 0x7888) {
        result.d4_orbits = orbit_sizes(d4_min);
        result.d5_orbits = orbit_sizes(d5_min);
    }
    return result;
}

void print_code(const Code& code) {
    char separator = '[';
    for (Mask v : code) { std::cout << separator << unsigned(v); separator = ','; }
    if (code.empty()) std::cout << '[';
    std::cout << ']';
}

void print_hist(const Hist& hist) {
    std::cout << '[' << hist[0];
    for (size_t i = 1; i < hist.size(); ++i) std::cout << ',' << hist[i];
    std::cout << ']';
}

void print_orbits(const std::map<int, int>& counts) {
    std::cout << '{';
    const char* sep = "";
    for (auto [size, count] : counts) {
        std::cout << sep << '"' << size << "\":" << count;
        sep = ",";
    }
    std::cout << '}';
}

int main(int argc, char** argv) {
    unsigned nthreads = argc == 2 ? unsigned(std::stoul(argv[1])) : std::thread::hardware_concurrency();
    if (!nthreads) return 2;
    const auto start = std::chrono::steady_clock::now();
    auto chars = characters();
    Code affine(chars.begin(), chars.end());
    for (Mask m : chars) affine.push_back(m ^ 0xffff);
    std::sort(affine.begin(), affine.end());
    Code qb;
    for (unsigned i = 0; i < 4; ++i) {
        for (unsigned j = i + 1; j < 4; ++j) qb.push_back(chars[1 << i] & chars[1 << j]);
    }
    auto quad = span(qb);
    Code seeds, hyperplanes;
    for (Mask q : quad) if (bent(q, chars)) seeds.push_back(q);
    for (Mask a : affine) if (std::popcount(a) == 8) hyperplanes.push_back(a);
    auto affine_spaces = subspaces(affine, 5);
    std::vector<std::vector<Kernel>> kernels(6);
    for (int t = 0; t <= 5; ++t) {
        for (const auto& space : affine_spaces[t]) {
            Kernel kernel;
            kernel.basis = basis(space);
            for (Mask w : space) kernel.support |= w;
            std::set<Mask> reps;
            for (Mask a : affine) {
                Mask least = a;
                for (Mask w : space) least = std::min(least, Mask(a ^ w));
                reps.insert(least);
            }
            kernel.coset_representatives.assign(reps.begin(), reps.end());
            kernels[t].push_back(kernel);
        }
    }
    std::vector<Result> results(seeds.size());
    std::atomic<size_t> job{0};
    std::vector<std::thread> workers;
    for (unsigned i = 0; i < nthreads; ++i) {
        workers.emplace_back([&] {
            while (true) {
                size_t j = job.fetch_add(1);
                if (j >= seeds.size()) return;
                results[j] = run_seed(seeds[j], quad, affine, chars, kernels, hyperplanes);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    std::cout << "{\"worker_threads\":" << nthreads << ",\"bent_quadratic_bases\":" << seeds.size()
              << ",\"rows\":[";
    for (size_t j = 0; j < seeds.size(); ++j) {
        const auto& r = results[j];
        if (j) std::cout << ',';
        if (r.census != results[0].census || r.with_e4 != results[0].with_e4) std::abort();
        std::cout << "{\"seed\":" << seeds[j] << ",\"legal\":" << r.legal << ",\"dimensions\":{";
        for (int d = 4; d <= 7; ++d) {
            if (d != 4) std::cout << ',';
            std::cout << '"' << d << "\":{\"hist\":";
            print_hist(r.census[d]);
            std::cout << ",\"with_e4\":";
            print_hist(r.with_e4[d]);
            std::cout << ",\"minimum_and_basis\":";
            print_code(r.witness[d]);
            std::cout << '}';
        }
        std::cout << "},\"d4_orbits\":";
        print_orbits(r.d4_orbits);
        std::cout << ",\"d5_support_le_15_orbits\":";
        print_orbits(r.d5_orbits);
        std::cout << '}';
    }
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::cout << "],\"seconds\":" << elapsed.count() << "}\n";
}
