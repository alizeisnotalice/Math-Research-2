// Exact marked simple-path counts in early first-contact words.
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

using Mask = std::uint32_t;
using Count = __uint128_t;
using Key = std::tuple<int,int,int,int>;
using Table = std::map<Key,Count>;
constexpr int MaxN=25;

Mask bit(int i) { return Mask{1}<<i; }
int first(Mask a) { return std::countr_zero(a); }
void drop(Mask& a) { a&=a-1; }
std::string decimal(Count value) {
    std::string out;
    do { out.push_back('0'+static_cast<int>(value%10)); value/=10; } while(value);
    std::reverse(out.begin(),out.end());
    return out;
}

class Census {
public:
    explicit Census(int n):n_(n),N_(n*n),all_(bit(N_)-1),ring_(all_+1),suffix_(all_+1) {
        for (int v=0;v<N_;++v) {
            for (int u=0;u<N_;++u) if (distance(v,u)==1) edges_[v]|=bit(u);
        }
        for (Mask s=1;s<=all_;++s) ring_[s]=ring_[s&(s-1)]|edges_[first(s)];
        suffix_[all_]=1;
        for (Mask s=all_-1;s>0;--s) {
            for (Mask z=ring_[s]&~s;z;drop(z)) {
                suffix_[s]+=suffix_[s|bit(first(z))];
            }
        }
    }

    std::vector<std::pair<int,int>> orbits() const {
        std::vector<std::pair<int,int>> out;
        for (int q=0;q<N_;++q) {
            std::set<int> images;
            for (int mode=0;mode<8;++mode) {
                int x=q%n_,y=q/n_;
                if (mode&1) std::swap(x,y);
                if (mode&2) x=n_-1-x;
                if (mode&4) y=n_-1-y;
                images.insert(y*n_+x);
            }
            if (q==*images.begin()) out.emplace_back(q,images.size());
        }
        return out;
    }

    Table one(int q,int orbit) const {
        std::vector<int> verts;
        Mask allowed=all_&~(bit(q)|edges_[q]);
        for (Mask a=allowed;a;drop(a)) verts.push_back(first(a));
        const Mask size=bit(static_cast<int>(verts.size()));
        std::vector<Mask> expanded(size),fringe(size);
        std::array<Mask,MaxN> adj{};
        for (int i=0;i<static_cast<int>(verts.size());++i) {
            for (int j=0;j<static_cast<int>(verts.size());++j) {
                if (edges_[verts[i]]&bit(verts[j])) adj[i]|=bit(j);
            }
        }
        for (Mask s=1;s<size;++s) {
            const int v=first(s);
            const Mask tail=s&(s-1);
            expanded[s]=expanded[tail]|bit(verts[v]);
            fringe[s]=fringe[tail]|adj[v];
        }
        std::vector<Count> weights(size), ways(size);
        Table out;
        const std::set<int> thresholds{std::max(3,n_/2+1),n_-1};
        for (int p=0;p<static_cast<int>(verts.size());++p) {
            const int P=verts[p], d=distance(P,q);
            if (d>=*thresholds.rbegin()) continue;
            std::fill(ways.begin(),ways.end(),0);
            ways[bit(p)]=1;
            for (Mask s=bit(p);s<size;++s) {
                if (!ways[s]) continue;
                const int len=std::popcount(s);
                if (len<6) {
                    for (Mask v=fringe[s]&~s;v;drop(v)) ways[s|bit(first(v))]+=ways[s];
                }
                const Mask K=expanded[s];
                const Mask U=all_&~(K|ring_[K]);
                if (!weights[s]) {
                    for (Mask v=edges_[q];v;drop(v)) weights[s]+=suffix_[K|bit(q)|bit(first(v))];
                    assert(weights[s]>0);
                }
                const Count coefficient=ways[s]*len*orbit*weights[s];
                for (int t:thresholds) {
                    if (d>=t) continue;
                    std::array<std::uint64_t,5> counts{};
                    Mask choices=U;
                    for (Mask a=U;a;drop(a)) {
                        const int v=first(a);
                        if (distance(P,v)<t && (distance(q,v)%2)==0) choices&=~bit(v);
                    }
                    paths(q,choices,P,t,0,counts);
                    for (int k=1;k<=4;++k) {
                        if (counts[k]) out[{len+1,d,t,k}]+=coefficient*counts[k];
                    }
                }
            }
        }
        return out;
    }

private:
    void paths(int at,Mask remain,int p,int threshold,int length,
               std::array<std::uint64_t,5>& counts) const {
        if (length && distance(p,at)>=threshold) ++counts[length];
        if (length==4) return;
        for (Mask v=edges_[at]&remain;v;drop(v)) {
            const int u=first(v);
            paths(u,remain&~bit(u),p,threshold,length+1,counts);
        }
    }

    int distance(int p,int q) const { return std::abs(p%n_-q%n_)+std::abs(p/n_-q/n_); }
    int n_,N_;
    Mask all_;
    std::vector<Mask> ring_;
    std::vector<Count> suffix_;
    std::array<Mask,MaxN> edges_{};
};

int main(int argc,char** argv) {
    assert(argc==3 || argc==5);
    const int n=std::stoi(argv[1]),threads=std::stoi(argv[2]);
    assert(2<=n && n<=5 && threads>=1);
    Census census(n);
    const auto tasks=census.orbits();
    const auto lo=argc==5 ? std::stoul(argv[3]) : 0;
    const auto hi=argc==5 ? std::stoul(argv[4]) : tasks.size();
    assert(lo<hi && hi<=tasks.size());
    const int workers=std::min(threads,static_cast<int>(hi-lo));
    std::atomic_size_t next{lo};
    std::vector<Table> partial(workers);
    std::vector<std::thread> pool;
    for (int i=0;i<workers;++i) {
        pool.emplace_back([&,i]{
            for (std::size_t j=next++;j<hi;j=next++) {
                const auto [q,orbit]=tasks[j];
                const auto table=census.one(q,orbit);
                for (const auto& [key,row]:table) {
                    partial[i][key]+=row;
                }
            }
        });
    }
    for (auto& thread:pool) thread.join();
    Table rows;
    for (const auto& table:partial) {
        for (const auto& [key,row]:table) {
            rows[key]+=row;
        }
    }
    for (const auto& [key,row]:rows) {
        const auto [rank,d,t,k]=key;
        std::cout<<rank<<' '<<d<<' '<<t<<' '<<k<<' '<<decimal(row)<<'\n';
    }
}
