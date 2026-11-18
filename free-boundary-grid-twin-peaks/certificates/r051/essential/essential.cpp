// Exact first essential contact after merging a low branch of size at most eight.
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

using Mask=std::uint32_t;
using Count=__uint128_t;
using Row=std::array<Count,9>;
using Cells=std::map<std::tuple<int,int,int,int,int>,Count>;
constexpr int MaxN=25;
Mask bit(int v) { return Mask{1}<<v; }
int first(Mask m) { return std::countr_zero(m); }
void drop(Mask& m) { m&=m-1; }
std::string decimal(Count x) {
    std::string out;
    do { out.push_back('0'+static_cast<int>(x%10)); x/=10; } while(x);
    std::reverse(out.begin(),out.end());
    return out;
}

struct Compact {
    std::vector<int> expanded_vertices;
    std::array<Mask,MaxN> compressed{};
    std::vector<Mask> adjacency;

    Compact(Mask vertices,const std::array<Mask,MaxN>& edges) {
        for (Mask v=vertices;v;drop(v)) {
            compressed[first(v)]=bit(static_cast<int>(expanded_vertices.size()));
            expanded_vertices.push_back(first(v));
        }
        for (int v:expanded_vertices) {
            adjacency.push_back(compress(edges[v]));
        }
    }
    Mask compress(Mask global) const {
        Mask out=0;
        for (Mask t=global;t;drop(t)) {
            out|=compressed[first(t)];
        }
        return out;
    }
    Mask size() const { return bit(static_cast<int>(expanded_vertices.size())); }
};

std::vector<Mask> neighborhoods(const Compact& region) {
    std::vector<Mask> out(region.size());
    for (Mask s=1;s<region.size();++s) {
        int v=first(s);
        out[s]=out[s^bit(v)]|region.adjacency[v];
    }
    return out;
}

std::vector<Count> completions(const Compact& region) {
    const Mask full=region.size()-1;
    const auto around=neighborhoods(region);
    std::vector<Count> ways(region.size());
    ways[full]=1;
    for (Mask state=full-1;state;--state) {
        for (Mask frontier=around[state]&~state;frontier;drop(frontier)) {
            ways[state]+=ways[state|bit(first(frontier))];
        }
    }
    return ways;
}

class Counter {
public:
    Counter(int n,int lo,int hi):n_(n),N_(n*n),min_branch_(lo),max_branch_(hi),full_(bit(N_)-1) {
        for (int v=0;v<N_;++v) {
            for (int u=0;u<N_;++u) {
                if (distance(u,v)==1) { edges_[v]|=bit(u); }
            }
            for (int sym=0;sym<8;++sym) {
                int x=v%n_,y=v/n_;
                if (sym&1) { std::swap(x,y); }
                if (sym&2) { x=n_-1-x; }
                if (sym&4) { y=n_-1-y; }
                symmetries_[sym][v]=y*n_+x;
            }
        }
    }

    void run(int workers,bool plan) {
        const auto tasks=work();
        assert(tasks.size()>0);
        std::uint64_t cap=(std::uint64_t{1}<<N_)*(sizeof(Count)+sizeof(Mask));
        std::vector<std::uint64_t> per;
        for (const auto& [branch,orbit]:tasks) {
            const int b=std::popcount(branch);
            const int free= N_-std::popcount(branch|boundary(branch));
            const int rim=std::popcount(boundary(branch));
            const auto bytes=(std::uint64_t{1}<<(N_-b))*(2*sizeof(Count)+2*sizeof(Mask))+
                             (std::uint64_t{1}<<free)*(b*sizeof(Row)+sizeof(Mask)+16*rim);
            per.push_back(bytes);
            std::cerr<<"task "<<branch<<" b "<<b<<" orbit "<<orbit<<" free "<<free<<" planned "<<bytes<<'\n';
        }
        std::sort(per.rbegin(),per.rend());
        for (int j=0;j<workers && j<static_cast<int>(per.size());++j) { cap+=per[j]; }
        std::cerr<<"aggregate array upper bytes "<<cap<<'\n';
        if (plan) { return; }
        assert(cap<(8ULL<<30));
        all_=completions(Compact(full_,edges_));
        std::atomic_size_t next{0};
        std::vector<Cells> results(workers);
        std::vector<std::thread> threads;
        for (int j=0;j<workers;++j) {
            threads.emplace_back([&,j]{
                for (std::size_t i=next++;i<tasks.size();i=next++) {
                    add_branch(tasks[i].first,tasks[i].second,results[j]);
                }
            });
        }
        for (auto& thread:threads) { thread.join(); }
        Cells result;
        for (const auto& local:results) {
            for (const auto& [key,x]:local) {
                result[key]+=x;
            }
        }
        for (const auto& [key,x]:result) {
            const auto [i,b,d,ell,h]=key;
            if (x) {
                std::cout << i << ' ' << b << ' ' << d << ' ' << ell << ' ' << h
                          << ' ' << decimal(x) << '\n';
            }
        }
    }
private:
    int distance(int u,int v) const { return std::abs(u%n_-v%n_)+std::abs(u/n_-v/n_); }
    Mask boundary(Mask S) const {
        Mask union_=0;
        for (Mask s=S;s;drop(s)) { union_|=edges_[first(s)]; }
        return union_&~S;
    }
    std::vector<std::pair<Mask,int>> work() const {
        std::vector<std::pair<Mask,int>> out;
        std::set<Mask> seen,layer;
        for (int v=0;v<N_;++v) { layer.insert(bit(v)); }
        for (int b=1;b<=max_branch_;++b) {
            std::set<Mask> following;
            for (Mask B:layer) {
                const Mask D=boundary(B);
                if (b<max_branch_) {
                    for (Mask s=D;s;drop(s)) { following.insert(B|bit(first(s))); }
                }
                if (b<min_branch_ || (B|D)==full_ || seen.contains(B)) { continue; }
                const auto start=seen.size();
                for (const auto& sym:symmetries_) {
                    Mask image=0;
                    for (Mask s=B;s;drop(s)) { image|=bit(sym[first(s)]); }
                    seen.insert(image);
                }
                out.emplace_back(B,static_cast<int>(seen.size()-start));
            }
            layer=std::move(following);
        }
        return out;
    }
    Count rooted(Mask S,int root) const {
        if (S==bit(root)) { return 1; }
        Count total=0;
        for (Mask rest=S^bit(root);rest;drop(rest)) {
            const auto v=first(rest);
            if (edges_[v]&S) { total+=rooted(S^bit(v),root); }
        }
        return total;
    }

    static int choose(int n,int k) {
        int total=1;
        for (int j=1;j<=k;++j) { total=total*(n-j+1)/j; }
        return total;
    }


    static Mask expand_without(Mask s,Mask B) {
        for (Mask rem=B;rem;drop(rem)) {
            const Mask low=bit(first(rem))-1;
            s=(s&low)|((s&~low)<<1);
        }
        return s;
    }

    void add_branch(Mask B,int orbit,Cells& out) const {
        const int b=std::popcount(B);
        const Mask D=boundary(B);
        const Compact after(full_&~B,edges_);
        const Compact core(full_&~(B|D),edges_);
        const Mask size=after.size();
        const Mask boundary_mask=after.compress(D);
        const Mask cap=core.size()-1;
        const auto adj=neighborhoods(after);
        std::vector<Mask> embed(core.size());
        std::vector<int> roots;
        for (Mask root=B;root;drop(root)) {
            roots.push_back(first(root));
        }
        for (Mask a=1;a<=cap;++a) {
            int v=first(a);
            embed[a]=embed[a^bit(v)]|after.compressed[core.expanded_vertices[v]];
        }

        std::array<int,26> shuffles{};
        std::array<Mask,MaxN> touches{};
        std::array<std::array<int,MaxN>,MaxN> distance_to_contact{};
        std::array<Count,MaxN> low_orders{};
        for (int q:roots) {
            low_orders[q]=rooted(B,q);
        }
        for (int a=1;a<=std::popcount(cap);++a) {
            shuffles[a]=choose(a+b-1,b);
        }
        for (Mask w=D;w;drop(w)) {
            int v=first(w);
            touches[v]=core.compress(edges_[v]);
            Mask reached=0;
            Mask layer=edges_[v]&B;
            for (int length=0;layer;++length) {
                for (Mask f=layer;f;drop(f)) {
                    distance_to_contact[v][first(f)]=length;
                }
                reached|=layer;
                layer=boundary(reached)&B;
            }
            assert(reached==B);
        }
        std::vector<std::vector<Row>> coefficients;
        for (int q:roots) {
            auto& prefix=coefficients.emplace_back(core.size());
            for (Mask a=1;a<=cap;++a) {
                if (std::has_single_bit(a)) {
                    prefix[a][distance(q,core.expanded_vertices[first(a)])]=1;
                    continue;
                }
                for (Mask choices=a;choices;drop(choices)) {
                    int v=first(choices);
                    if (core.adjacency[v]&a) {
                        for (int d=2;d<=2*n_-2;++d) {
                            prefix[a][d]+=prefix[a^bit(v)][d];
                        }
                    }
                }
            }
        }

        using Small=std::array<std::array<std::array<Count,8>,9>,26>;
        Small total{}, assigned{};
        std::vector<std::array<Mask,2>> initial;
        for (Mask a=1;a<=cap;++a) {
            bool has_words=false;
            for (int d=2;d<=2*n_-2;++d) {
                has_words|=coefficients[0][a][d]>0;
            }
            if (!has_words) { continue; }
            for (Mask ws=D;ws;drop(ws)) {
                int w=first(ws);
                if (a&touches[w]) {
                    initial.push_back({a,static_cast<Mask>(w)});
                }
            }
        }
        auto accumulate=[&](Small& data,Mask a,int w,Count suffix) {
            if (!suffix) { return; }
            int cardinality=std::popcount(a);
            int i=cardinality+b+1;
            Count factor=Count(orbit)*shuffles[cardinality]*suffix;
            for (std::size_t k=0;k<roots.size();++k) {
                int q=roots[k];
                const int ell=distance_to_contact[w][q];
                for (int d=2;d<=2*n_-2;++d) {
                    data[i][d][ell]+=coefficients[k][a][d]*factor*low_orders[q];
                }
            }
        };
        auto initial_state=[&](Mask a,int w) {
            return embed[a]|after.compressed[w];
        };
        for (auto [a,w]:initial) {
            Mask s=initial_state(a,w);
            accumulate(total,a,w,all_[expand_without(s,B)|B]);
        }

        std::vector<Count> prev(size),current(size);
        std::vector<Mask> live;
        for (Mask s=1;s<size;++s) {
            Mask essential=boundary_mask&~(s|adj[s]);
            if (!essential) { continue; }
            live.push_back(s);
            for (Mask e=essential;e;drop(e)) {
                prev[s]+=all_[expand_without(s|(e&-e),B)|B];
            }
        }
        for (int h=1;h<=N_-b-2;++h) {
            Small layer{};
            for (auto [a,w]:initial) {
                Mask s=initial_state(a,w);
                accumulate(layer,a,w,prev[s]);
            }
            for (int i=b+2;i<=N_;++i) {
                for (int d=2;d<=2*n_-2;++d) {
                    for (int ell=0;ell<=b-1;++ell) {
                        Count value=layer[i][d][ell];
                        if (value) {
                            assert(h<=N_-i);
                            out[{i,b,d,ell,h}]+=value;
                            assigned[i][d][ell]+=value;
                        }
                    }
                }
            }
            bool next_nonzero=false;
            for (Mask s:live) {
                Count value=0;
                if (std::popcount(s)+h < N_-b) {
                    for (Mask f=adj[s]&~s;f;drop(f)) {
                        value+=prev[s|(f&-f)];
                    }
                }
                current[s]=value;
                next_nonzero|=value!=0;
            }
            if (!next_nonzero) { break; }
            prev.swap(current);
        }
        for (int i=b+2;i<=N_;++i) {
            for (int d=2;d<=2*n_-2;++d) {
                for (int ell=0;ell<=b-1;++ell) {
                    const Count a=total[i][d][ell],c=assigned[i][d][ell];
                    assert(c<=a);
                    if (a>c) {
                        out[{i,b,d,ell,N_+1-i}]+=a-c;
                    }
                }
            }
        }
    }

    const int n_,N_,min_branch_,max_branch_;
    const Mask full_;
    std::array<Mask,MaxN> edges_{};
    std::array<std::array<int,MaxN>,8> symmetries_{};
    std::vector<Count> all_;
};

int main(int argc,char** argv) {
    assert(argc==5||argc==6);
    const int n=std::atoi(argv[1]),workers=std::atoi(argv[2]);
    const int lo=std::atoi(argv[3]),hi=std::atoi(argv[4]);
    assert(n>=3 && n<=5 && workers>=1 && workers<=12 && lo>=1 && lo<=hi && hi<=8);
    const bool plan=argc==6 && std::string(argv[5])=="--plan";
    assert(argc==5||plan);
    Counter(n,lo,hi).run(workers,plan);
}
