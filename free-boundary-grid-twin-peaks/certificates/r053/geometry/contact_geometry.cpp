// Exact census of the vertices reachable before the low peak's first neighbor.
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
using Key = std::tuple<int,int,int>;
using Row = std::array<Count,2>;
using Table = std::map<Key,Row>;
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
        std::vector<Row> weights(size);
        // Four extrema pack into four nibbles on the supported grids.
        std::vector<std::uint16_t> extents(size);
        std::vector<Count> ways(size);
        Table out;
        for (int p=0;p<static_cast<int>(verts.size());++p) {
            std::fill(ways.begin(),ways.end(),0);
            const int P=verts[p],d=distance(P,q);
            ways[bit(p)]=1;
            for (Mask s=bit(p);s<size;++s) {
                if (!ways[s]) continue;
                for (Mask v=fringe[s]&~s;v;drop(v)) {
                    ways[s|bit(first(v))]+=ways[s];
                }
                if (!weights[s][0]) {
                    const Mask K=expanded[s];
                    const Mask U=all_&~(K|ring_[K]);
                    for (Mask v=edges_[q];v;drop(v)) {
                        const Mask child=bit(first(v));
                        const Count count=suffix_[K|bit(q)|child];
                        weights[s][0]+=count;
                        if (U&child) weights[s][1]+=count;
                    }
                    assert(weights[s][0]>0);
                    Mask comp=bit(q);
                    while (true) {
                        const Mask next=comp|(ring_[comp]&U);
                        if (next==comp) break;
                        comp=next;
                    }
                    int lo1=15,hi1=0,lo2=15,hi2=0;
                    for (Mask a=comp;a;drop(a)) {
                        const int x=first(a)%n_,y=first(a)/n_;
                        const int u=x+y,v=x-y+n_-1;
                        lo1=std::min(lo1,u); hi1=std::max(hi1,u);
                        lo2=std::min(lo2,v); hi2=std::max(hi2,v);
                    }
                    extents[s]=static_cast<std::uint16_t>(lo1|(hi1<<4)|(lo2<<8)|(hi2<<12));
                }
                const int extent=extents[s],u=P%n_+P/n_,v=P%n_-P/n_+n_-1;
                const int far=std::max({u-(extent&15),((extent>>4)&15)-u,
                                    v-((extent>>8)&15),(extent>>12)-v});
                assert(far>=d);
                const int len=std::popcount(s);
                Row& row=out[{len+1,d,far}];
                const Count coefficient=ways[s]*len*orbit;
                row[0]+=coefficient*weights[s][0];
                row[1]+=coefficient*weights[s][1];
            }
        }
        return out;
    }

private:
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
                    for (int a=0;a<2;++a) partial[i][key][a]+=row[a];
                }
            }
        });
    }
    for (auto& thread:pool) thread.join();
    Table rows;
    for (const auto& table:partial) {
        for (const auto& [key,row]:table) {
            for (int a=0;a<2;++a) rows[key][a]+=row[a];
        }
    }
    for (const auto& [key,row]:rows) {
        const auto [rank,d,far]=key;
        std::cout<<rank<<' '<<d<<' '<<far<<' '<<decimal(row[0])<<' '<<decimal(row[1])<<'\n';
    }
}
