// Independent exact C-budget prefix traversal for the S={3} model.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using I64 = std::int64_t;
using I128 = __int128_t;
using Row = std::array<I64,9>;
constexpr I64 cofactor_cap = I64{1} << 40;
constexpr I64 convolution_cap = I64{1} << 56;
constexpr int depth_cap = 96;
constexpr std::size_t state_cap = 10000;

struct Ring {
    I128 u, v; // (u + v sqrt(D))/2.
};

struct Interval {
    I128 min, max;
};

I128 floor_half(I128 v) {
    return v/2 - (v<0 && v%2!=0);
}

I128 sqrt_floor(I128 n) {
    // For these bounded operands even the guessed square stays inside I128.
    I128 k=static_cast<I64>(std::sqrt(static_cast<double>(n)));
    for (; k*k>n; --k) {}
    for (; (k+1)*(k+1)<=n; ++k) {}
    return k;
}

Interval allowed(Ring b, I64 D, int row, bool reversed) {
    const Interval empty{1,0};
    if ((!reversed && row==5) || (reversed && row==3)) {
        if (b.v!=0) return empty;
        if (b.u%2!=0) throw std::runtime_error("nonintegral rational base");
        return {-b.u/2,-b.u/2};
    }
    if (!reversed && row==2) {
        if (b.v!=0) throw std::runtime_error("unexpected irrational row 2");
        return b.u==0 || b.u==2 ? Interval{0,0} : empty;
    }
    if (b.v==0) {
        if (b.u%2!=0) throw std::runtime_error("nonintegral rational base");
        return {-b.u/2,1-b.u/2};
    }
    I128 r=sqrt_floor(b.v*b.v*D);
    if (b.v<0) r=-r-1;
    const I128 g=-floor_half(b.u+r);
    return {g,g};
}

struct State {
    std::array<I64,depth_cap+1> G{}, Q{}, W{};
    State() { G[0]=Q[0]=W[0]=1; }
};

Ring partial_convolution(const std::array<Ring,9>& p, const std::array<I64,depth_cap+1>& g, int k) {
    Ring u{};
    for (int j=1;j<=std::min(8,k);++j) {
        u.u+=p[j].u*g[k-j];
        u.v+=p[j].v*g[k-j];
    }
    return u;
}

struct Single {
    bool empty;
    int layers;
    std::size_t max_live;
    const char* reason;
};

Single traverse(const Row& row, bool reversed) {
    const auto [L,M,T,N,D,b,c,d,e]=row;
    const int pos=reversed ? 5 : 3;
    std::array<I64,9> f{1,L,M,0,b,c,d,e,1};
    if (reversed) std::reverse(f.begin(),f.end());
    std::array<Ring,9> p{}, pb{};
    std::array<I64,17> H{};
    for (int j=0;j<=8;++j) {
        for (int k=0;k<=8;++k) H[j+k]+=f[j]*f[k];
        H[j+pos]+=T*f[j];
        p[j]={2*f[j],0};
        pb[j]=p[j];
    }
    H[2*pos]+=N;
    p[pos]={T,1};
    pb[pos]={T,-1};
    std::vector<State> states, next;
    states.reserve(32);
    next.reserve(32);
    states.emplace_back();
    std::size_t max_live=0;
    for (int k=1;k<=depth_cap;++k) {
        next.clear();
        for (const State& s : states) {
            const Interval gl=allowed(partial_convolution(p,s.G,k),D,k,reversed);
            const Interval ql=allowed(partial_convolution(pb,s.Q,k),D,k,reversed);
            if (gl.min>gl.max || ql.min>ql.max) continue;
            I128 internal=0;
            for (int j=1;j<k;++j)
                internal+=I128(s.G[j])*s.Q[k-j];
            I128 cross=0;
            for (int j=1;j<=std::min(16,k);++j)
                cross+=I128(H[j])*s.W[k-j];
            for (int C=0;C<=1;++C) {
                const I128 sum_gq=C-cross-internal;
                const I128 lo=std::max(gl.min,sum_gq-ql.max);
                const I128 hi=std::min(gl.max,sum_gq-ql.min);
                for (I128 g=lo;g<=hi;++g) {
                    const I128 q=sum_gq-g;
                    const I128 w=C-cross;
                    if (g<-cofactor_cap || g>cofactor_cap || q<-cofactor_cap ||
                        q>cofactor_cap || w<-convolution_cap || w>convolution_cap)
                        return {false,k,max_live,"arithmetic_cap"};
                    next.push_back(s);
                    next.back().G[k]=static_cast<I64>(g);
                    next.back().Q[k]=static_cast<I64>(q);
                    next.back().W[k]=static_cast<I64>(w);
                    if (next.size()>state_cap) return {false,k,max_live,"state_cap"};
                }
            }
        }
        if (next.empty()) return {true,k,max_live,"empty"};
        max_live=std::max(max_live,next.size());
        states.swap(next);
    }
    return {false,depth_cap,max_live,"depth_cap"};
}

Row read_row(std::string line) {
    for (auto& ch : line) if (ch=='[' || ch==']' || ch==',') ch=' ';
    std::istringstream input(std::move(line));
    Row row{};
    for (auto& item : row) if (!(input>>item)) throw std::runtime_error("bad input row");
    std::string tail;
    if (input>>tail) throw std::runtime_error("extra input column");
    const auto [L,M,T,N,D,b,c,d,e]=row;
    if (L<-4 || L>4 || M!=L*L || T<-104 || T>104 || N<-117 || N>2702
        || D<5 || D>468 || T*T-4*N!=D || sqrt_floor(D)*sqrt_floor(D)==D
        || b<-120 || b>120 || c<-92 || c>92 || d<-41 || d>41 || e<-9 || e>9)
        throw std::runtime_error("input row out of bounds");
    return row;
}

int main(int argc, char** argv) {
    if (argc!=2) throw std::runtime_error("usage: s3_budget_tree ROOT_PASSES_JSONL");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open input");
    std::array<std::uint64_t,depth_cap+1> death_f{},death_r{};
    std::uint64_t count=0, unknowns=0, all_layers=0;
    std::size_t biggest=0;
    for (std::string line; std::getline(input,line);) {
        const Row row=read_row(line);
        const auto forward=traverse(row,false);
        const auto reverse=forward.empty ? Single{false,0,0,"not_needed"} : traverse(row,true);
        all_layers+=forward.layers+reverse.layers;
        biggest=std::max({biggest,forward.max_live,reverse.max_live});
        ++count;
        const char status=forward.empty ? 'F' : (reverse.empty ? 'R' : '?');
        if(status=='F') ++death_f[forward.layers];
        else if(status=='R') ++death_r[reverse.layers];
        else ++unknowns;
        std::cout << '[' << row[0];
        for (std::size_t i=1;i<row.size();++i) std::cout << ',' << row[i];
        std::cout << "],\"" << status << "\"," << forward.layers << ',' << reverse.layers << ','
                  << (reverse.empty ? reverse.max_live : forward.max_live);
        if (status=='?') std::cout << ",\"" << forward.reason << "\",\"" << reverse.reason << '"';
        std::cout << '\n';
    }
    if (!input.eof()) throw std::runtime_error("input read error");
    std::cerr << "{\"count\":" << count << ",\"layers\":" << all_layers << ",\"unknowns\":" << unknowns
              << ",\"largest\":" << biggest << ",\"forward\":{";
    const auto print=[&](const auto& values) {
        bool first=true;
        for (int k=1;k<=depth_cap;++k) if (values[k]) {
            std::cerr << (first ? "" : ",") << '"' << k << "\":" << values[k];
            first=false;
        }
    };
    print(death_f);
    std::cerr << "},\"reverse\":{";
    print(death_r);
    std::cerr << "}}\n";
}
