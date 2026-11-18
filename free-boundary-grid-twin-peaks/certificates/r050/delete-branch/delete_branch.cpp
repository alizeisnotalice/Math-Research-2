// Count original words for which deleting the fixed low branch keeps one peak.
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
using Cells=std::map<std::tuple<int,int,int>,std::array<Count,2>>;
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
    Counter(int n,int cap):n_(n),N_(n*n),max_branch_(cap),full_(bit(N_)-1) {
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
            const auto bytes=(std::uint64_t{1}<<(N_-b))*(sizeof(Count)+sizeof(Mask))+
                             (std::uint64_t{1}<<free)*(sizeof(Row)+sizeof(Mask));
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
                result[key][0]+=x[0];
                result[key][1]+=x[1];
            }
        }
        for (const auto& [key,x]:result) {
            const auto [i,b,d]=key;
            assert(x[1]<=x[0]);
            std::cout<<i<<' '<<b<<' '<<d<<' '<<decimal(x[0])<<' '<<decimal(x[1])<<'\n';
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
                if ((B|D)==full_||seen.contains(B)) { continue; }
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
        for (Mask avail=S^bit(root);avail;drop(avail)) {
            const int v=first(avail);
            if (edges_[v]&S) { total+=rooted(S^bit(v),root); }
        }
        return total;
    }
    static int choose(int n,int k) {
        int total=1;
        for (int j=1;j<=k;++j) { total=total*(n-j+1)/j; }
        return total;
    }

    void add_branch(Mask B,int orbit,Cells& out) const {
        const int b=std::popcount(B);
        const Mask D=boundary(B);
        const Compact deleted(full_&~B,edges_);
        const auto continuation=completions(deleted);
        const Compact source(full_&~(B|D),edges_);
        const Mask cap=source.size()-1;
        std::vector<Mask> expanded(source.size());
        for (Mask s=1;s<=cap;++s) {
            const int v=first(s);
            expanded[s]=expanded[s^bit(v)]|bit(source.expanded_vertices[v]);
        }
        std::array<int,26> shuffles{};
        for (int a=1;a<=std::popcount(cap);++a) {
            shuffles[a]=choose(a+b-1,b);
        }

        for (Mask roots=B;roots;drop(roots)) {
            const int q=first(roots);
            const Count low=rooted(B,q);
            assert(low>0);
            std::vector<Row> prefix(source.size());
            for (Mask s=1;s<=cap;++s) {
                if (std::has_single_bit(s)) {
                    const int p=source.expanded_vertices[first(s)];
                    prefix[s][distance(p,q)]=1;
                    continue;
                }
                for (Mask v=s;v;drop(v)) {
                    const int j=first(v);
                    if (source.adjacency[j]&s) {
                        for (int d=2;d<=2*n_-2;++d) {
                            prefix[s][d]+=prefix[s^bit(j)][d];
                        }
                    }
                }
            }
            for (Mask ws=D;ws;drop(ws)) {
                const int w=first(ws);
                const Mask attach=source.compress(edges_[w]);
                for (Mask s=1;s<=cap;++s) {
                    if (!(s&attach)) { continue; }
                    const int a=std::popcount(s),i=a+b+1;
                    const Mask kept=expanded[s]|bit(w);
                    const Count all=all_[kept|B];
                    const Count legal=continuation[deleted.compress(kept)];
                    assert(legal<=all);
                    const Count factor=orbit*low*shuffles[a];
                    for (int d=2;d<=2*n_-2;++d) {
                        const Count root_factor=prefix[s][d]*factor;
                        if (root_factor) {
                            auto& row=out[{i,b,d}];
                            row[0]+=root_factor*all;
                            row[1]+=root_factor*legal;
                        }
                    }
                }
            }
        }
    }

    const int n_,N_,max_branch_;
    const Mask full_;
    std::array<Mask,MaxN> edges_{};
    std::array<std::array<int,MaxN>,8> symmetries_{};
    std::vector<Count> all_;
};

int main(int argc,char** argv) {
    assert(argc==4||argc==5);
    const int n=std::atoi(argv[1]),workers=std::atoi(argv[2]),branch=std::atoi(argv[3]);
    assert(n>=3&&n<=5&&workers>=1&&workers<=12&&branch>=1&&branch<=8);
    const bool plan=argc==5&&std::string(argv[4])=="--plan";
    assert(argc==4||plan);
    Counter(n,branch).run(workers,plan);
}
