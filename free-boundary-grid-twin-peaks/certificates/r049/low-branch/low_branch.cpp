// Exact next-contact table for low branches with at most eight vertices.
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

using Mask = std::uint32_t;
using Count = __uint128_t;
using Row = std::array<Count, 9>;
using Key = std::tuple<int, int, int, int, int>;
using Histogram = std::map<Key, Count>;

Mask bit(int v) { return Mask{1} << v; }
int first(Mask s) { return std::countr_zero(s); }
void drop(Mask& s) { s &= s - 1; }

std::string decimal(Count x) {
    std::string s;
    do {
        s.push_back('0' + static_cast<int>(x % 10));
        x /= 10;
    } while (x);
    std::reverse(s.begin(), s.end());
    return s;
}

struct Task {
    Mask branch;
    int root, merge_point, multiplicity, available;
};

class Counter {
  public:
    Counter(int n, int max_branch) : n_(n), N_(n*n), full_(bit(N_)-1), max_branch_(max_branch) {
        for (int v = 0; v < N_; ++v) {
            for (int u = 0; u < N_; ++u) {
                if (distance(u, v) == 1) {
                    adj_[v] |= bit(u);
                }
            }
            for (int g = 0; g < 8; ++g) {
                int x = v%n_, y = v/n_;
                if (g&1) { std::swap(x, y); }
                if (g&2) { x = n_-1-x; }
                if (g&4) { y = n_-1-y; }
                sym_[g][v] = y*n_+x;
            }
        }
    }

    void run(int workers, bool plan) {
        auto work = tasks();
        std::sort(work.begin(), work.end(), [](const Task& a, const Task& b) {
            return a.available > b.available;
        });
        auto bound = (std::uint64_t(full_)+1)*(sizeof(Mask)+sizeof(Count));
        for (std::size_t j = 0; j < work.size(); ++j) {
            const auto [B,q,w,m,u] = work[j];
            const auto bytes = (std::uint64_t{1}<<u)*(2*sizeof(Mask)+sizeof(Row)+2*sizeof(Count));
            std::cerr << "task " << B << ' ' << q << ' ' << w << " orbit " << m
                      << " free " << u << " arrays " << bytes << '\n';
            if (j < static_cast<std::size_t>(workers)) {
                bound += bytes;
            }
        }
        std::cerr << "total array upper bytes " << bound << '\n';
        if (plan) {
            return;
        }
        assert(bound < 9ULL*(1ULL<<30));
        build_suffix();
        std::atomic_size_t next{0};
        std::vector<Histogram> rows(workers);
        std::vector<std::thread> pool;
        for (int j = 0; j < workers; ++j) {
            pool.emplace_back([&,j] {
                for (std::size_t i = next++; i < work.size(); i = next++) {
                    for (const auto& [key,c] : count(work[i])) {
                        rows[j][key] += c;
                    }
                }
            });
        }
        for (auto& th : pool) {
            th.join();
        }
        Histogram total;
        for (const auto& result : rows) {
            for (const auto& [key,c] : result) {
                total[key] += c;
            }
        }
        for (const auto& [key,c] : total) {
            const auto [i,b,h,d,T] = key;
            std::cout << i << ' ' << b << ' ' << h << ' ' << d << ' ' << T << ' ' << decimal(c) << '\n';
        }
    }

  private:
    int distance(int u, int v) const {
        return std::abs(u%n_-v%n_)+std::abs(u/n_-v/n_);
    }

    Mask boundary(Mask B) const {
        Mask neighbors = 0;
        for (Mask s = B; s; drop(s)) {
            neighbors |= adj_[first(s)];
        }
        return neighbors & ~B;
    }

    std::vector<Task> tasks() const {
        std::vector<Task> work;
        std::set<std::tuple<Mask,int,int>> seen;
        for (int q = 0; q < N_; ++q) {
            std::set<Mask> branches{bit(q)};
            std::set<Mask> layer{bit(q)};
            for (int b=1; b<max_branch_; ++b) {
                std::set<Mask> following;
                for (Mask B : layer) {
                    for (Mask f=boundary(B); f; drop(f)) {
                        following.insert(B|bit(first(f)));
                    }
                }
                branches.insert(following.begin(),following.end());
                layer=std::move(following);
            }
            for (Mask B : branches) {
                const Mask D = boundary(B);
                if ((B|D)==full_) { continue; }
                for (Mask ws = D; ws; drop(ws)) {
                    const int w = first(ws);
                    if (seen.contains({B,q,w})) {
                        continue;
                    }
                    const auto old = seen.size();
                    for (const auto& tr : sym_) {
                        Mask image = 0;
                        for (Mask s = B; s; drop(s)) {
                            image |= bit(tr[first(s)]);
                        }
                        seen.insert({image,tr[q],tr[w]});
                    }
                    work.push_back({B,q,w,static_cast<int>(seen.size()-old),N_-std::popcount(B|D)});
                }
            }
        }
        return work;
    }

    int branch_words(Mask B, int root) const {
        std::vector<int> vertices{root};
        for (Mask choices=B^bit(root); choices; drop(choices)) {
            vertices.push_back(first(choices));
        }
        const int n=static_cast<int>(vertices.size());
        std::array<Mask,8> edges{};
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                if (adj_[vertices[i]]&bit(vertices[j])) {
                    edges[i]|=bit(j);
                }
            }
        }
        std::array<int,1<<8> count{};
        count[1]=1;
        for (Mask s=3; s<bit(n); s+=2) {
            for (Mask choices=s^1; choices; drop(choices)) {
                const auto v=first(choices);
                if (edges[v]&s) {
                    count[s]+=count[s^bit(v)];
                }
            }
        }
        return count[bit(n)-1];
    }

    static int choose(int n, int k) {
        int result=1;
        for (int v=1; v<=k; ++v) {
            result=result*(n-v+1)/v;
        }
        return result;
    }

    void build_suffix() {
        around_.resize(full_+1);
        suffix_.resize(full_+1);
        for (Mask s = 1; s <= full_; ++s) {
            const int v = first(s);
            around_[s] = around_[s^bit(v)] | adj_[v];
        }
        suffix_[full_] = 1;
        for (Mask s = full_-1; s > 0; --s) {
            Count total = 0;
            for (Mask f = around_[s]&~s; f; drop(f)) {
                total += suffix_[s|bit(first(f))];
            }
            suffix_[s] = total;
        }
    }

    Histogram count(const Task& task) const {
        const auto [B,q,w,mult,available] = task;
        const int b = std::popcount(B);
        const Mask D = boundary(B);
        const Mask U = full_&~(B|D);
        const Mask fixed = B|bit(w);
        const Mask target = D^bit(w);
        const int target_size=std::popcount(target);
        assert(target);  // Every permitted B has at least two outside neighbors.
        std::vector<int> of;
        std::array<Mask,25> into{};
        for (Mask s = U; s; drop(s)) {
            const int v = first(s);
            into[v] = bit(static_cast<int>(of.size()));
            of.push_back(v);
        }
        assert(static_cast<int>(of.size()) == available);
        const Mask cap = bit(available)-1;
        std::array<Mask,25> compact_adj{};
        for (int v = 0; v < N_; ++v) {
            for (Mask s = adj_[v]&U; s; drop(s)) {
                compact_adj[v] |= into[first(s)];
            }
        }
        std::vector<Mask> embedded(1+cap), near(1+cap);
        for (Mask s = 1; s <= cap; ++s) {
            const int u = of[first(s)];
            const Mask rest = s^bit(first(s));
            embedded[s] = embedded[rest]|bit(u);
            near[s] = near[rest]|compact_adj[u];
        }
        Mask fixed_near = 0;
        for (Mask s = fixed; s; drop(s)) {
            fixed_near |= compact_adj[first(s)];
        }

        const int max_d = 2*n_-2;
        std::vector<Row> prefix(1+cap);
        for (Mask s = 1; s <= cap; ++s) {
            if (std::has_single_bit(s)) {
                prefix[s][distance(q,of[first(s)])] = 1;
                continue;
            }
            for (Mask f = s; f; drop(f)) {
                const auto j = first(f);
                if (compact_adj[of[j]]&s) {
                    const auto& source = prefix[s^bit(j)];
                    for (int d = 2; d <= max_d; ++d) {
                        prefix[s][d] += source[d];
                    }
                }
            }
        }

        std::vector<Count> layers[2];
        layers[0].resize(cap+1);
        layers[1].resize(cap+1);
        for (Mask s = 0; s <= cap; ++s) {
            const Mask state = embedded[s]|fixed;
            for (Mask f = around_[state]&target; f; drop(f)) {
                layers[0][s] += suffix_[state|bit(first(f))];
            }
        }
        Histogram output;
        const int low_words=branch_words(B,q);
        std::array<int,26> shuffles{};
        for (int a=1; a<=available; ++a) {
            shuffles[a]=choose(a+b-1,b);
        }
        const Mask joinable = compact_adj[w];
        for (int h = 1; h <= available+1; ++h) {
            const auto& prev = layers[(h-1)&1];
            auto& current = layers[h&1];
            for (Mask s = 1; s <= cap; ++s) {
                if (prev[s] && (s&joinable)) {
                    const int a = std::popcount(s);
                    const int i = a+b+1;
                    const Count factor = Count(mult)*shuffles[a]*low_words*prev[s];
                    for (int d = 2; d <= max_d; ++d) {
                        if (prefix[s][d]) {
                            output[{i,b,h,d,target_size}] += factor*prefix[s][d];
                        }
                    }
                }
            }
            if (h <= available) {
                for (Mask s = 0; s <= cap; ++s) {
                    Count total = 0;
                    if (std::popcount(s) <= available-h) {
                        for (Mask f = (near[s]|fixed_near)&~s; f; drop(f)) {
                            total += prev[s|bit(first(f))];
                        }
                    }
                    current[s] = total;
                }
            }
        }
        return output;
    }

    const int n_,N_;
    const Mask full_;
    const int max_branch_;
    std::array<Mask,25> adj_{};
    std::array<std::array<int,25>,8> sym_{};
    std::vector<Mask> around_;
    std::vector<Count> suffix_;
};

int main(int argc, char** argv) {
    assert(argc == 4 || argc == 5);
    int n = std::atoi(argv[1]), workers = std::atoi(argv[2]);
    assert(n >= 2 && n <= 5 && workers >= 1 && workers <= 12);
    const int max_branch=std::atoi(argv[3]);
    assert(max_branch>=1 && max_branch<=8);
    const bool plan=argc==5 && std::string(argv[4])=="--plan";
    assert(argc==4 || plan);
    Counter(n,max_branch).run(workers,plan);
}
