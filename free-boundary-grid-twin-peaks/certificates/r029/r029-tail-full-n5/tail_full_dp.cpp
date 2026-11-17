// Exact punctured stable-fibre L-distribution on P_n square P_n.
//
// Fixed hole v.  First build f(S), the number of legal completions from S:
// a new point x may be born only if x is a neighbour of the hole; otherwise
// it needs an earlier neighbour in S.  Then, for every exterior root p,
// enumerate the prefix before the first p-neighbour with a compressed forward
// DP.  All arithmetic is checked unsigned __int128.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u128 = unsigned __int128;
using u32 = std::uint32_t;
static constexpr u128 UMAX = ~u128(0);

static std::string dec(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
    std::reverse(s.begin(), s.end());
    return s;
}
static u128 addc(u128 a, u128 b) {
    if (UMAX - a < b) throw std::overflow_error("u128 addition");
    return a + b;
}
static u128 mulc(u128 a, u128 b) {
    if (a && b > UMAX / a) throw std::overflow_error("u128 multiplication");
    return a * b;
}

struct Grid {
    int n, N, M, hole;
    std::vector<int> orig_of, local_of;
    std::vector<u32> adj;
    u32 nv = 0, full = 0;
    Grid(int side, int v): n(side), N(side*side), M(N-1), hole(v), local_of(N,-1) {
        if (n < 2 || M > 31 || v < 0 || v >= N) throw std::runtime_error("bad grid/hole");
        for (int q=0; q<N; ++q) if (q != hole) {
            local_of[q] = int(orig_of.size()); orig_of.push_back(q);
        }
        adj.assign(M,0); full=(u32(1)<<M)-1;
        const int dx[4]={1,-1,0,0}, dy[4]={0,0,1,-1};
        for (int i=0;i<M;++i) {
            int q=orig_of[i], x=q%n, y=q/n;
            for(int k=0;k<4;++k) {
                int xx=x+dx[k], yy=y+dy[k];
                if(xx<0||xx>=n||yy<0||yy>=n) continue;
                int r=yy*n+xx;
                if(r==hole) { nv |= (u32(1)<<i); }
                else adj[i] |= (u32(1)<<local_of[r]);
            }
        }
    }
};

struct HoleResult {
    int v; std::vector<u128> hist; u128 B=0, A=0;
    std::vector<u128> rootB;
};

static HoleResult compute_hole(int n, int v) {
    Grid g(n,v);
    const std::size_t states=std::size_t(1) << g.M;
    // 4 bytes/state for the union of neighbours; 16 bytes/state for f.
    std::vector<u32> nb(states); nb[0]=0;
    for (u32 S=1; S<=g.full; ++S) {
        u32 b=S&(~S+1); // lowest set bit, valid in unsigned arithmetic
        int i=__builtin_ctz(b);
        nb[S]=nb[S^b] | g.adj[i];
    }
    std::vector<u128> f(states); f[g.full]=1;
    for (u32 S=g.full; S-- > 0;) {
        u32 avail=(g.nv | nb[S]) & (g.full ^ S);
        u128 z=0;
        while(avail) { u32 b=avail&(~avail+1); avail^=b; z=addc(z,f[S|b]); }
        f[S]=z;
    }

    HoleResult out; out.v=v; out.hist.assign(g.M+1,0); out.rootB.assign(g.M,0);
    const u32 exterior=g.full & ~g.nv;
    for (u32 qq=exterior; qq;) {
        u32 bp=qq&(~qq+1); qq^=bp; int p=__builtin_ctz(bp);
        const u32 allowed=g.full & ~(bp | g.adj[p]);
        std::vector<int> verts;
        for(u32 a=allowed;a;) { u32 b=a&(~a+1); a^=b; verts.push_back(__builtin_ctz(b)); }
        const int K=int(verts.size());
        std::vector<u128> ways(std::size_t(1)<<K); ways[0]=1;
        std::vector<u32> expand(std::size_t(1)<<K); expand[0]=bp;
        for (u32 I=1; I<(u32(1)<<K); ++I) {
            u32 b=I&(~I+1); expand[I]=expand[I^b] | (u32(1)<<verts[__builtin_ctz(b)]);
        }
        std::vector<u128> hp(g.M+1,0);
        for (u32 I=0; I<(u32(1)<<K); ++I) {
            if (!ways[I]) continue;
            u32 S=expand[I];
            u128 suffix=0;
            for(u32 a=g.adj[p];a;) { u32 b=a&(~a+1); a^=b; suffix=addc(suffix,f[S|b]); }
            int t=__builtin_popcount(S);
            hp[t]=addc(hp[t],mulc(ways[I],suffix));
            u32 eligible=(g.nv | nb[S]) & allowed & ~S;
            while(eligible) {
                u32 b=eligible&(~eligible+1); eligible^=b;
                int x=__builtin_ctz(b);
                // locate x in the compressed coordinate set; K<=21, so a
                // tiny scan is preferable to another 2^K table.
                int j=0; while (verts[j]!=x) ++j;
                ways[I | (u32(1)<<j)]=addc(ways[I | (u32(1)<<j)],ways[I]);
            }
        }
        u128 bpB=0;
        for (int t=1;t<=g.M;++t) { out.hist[t]=addc(out.hist[t],hp[t]); bpB=addc(bpB,hp[t]); }
        if (bpB != f[bp]) throw std::runtime_error("prefix/suffix partition failed");
        out.rootB[p]=bpB;
    }
    for (int t=1;t<=g.M;++t) { out.B=addc(out.B,out.hist[t]); out.A=addc(out.A,mulc(t,out.hist[t])); }
    return out;
}

static void emit(const HoleResult& r, int n, bool comma) {
    int v=r.v;
    std::cout << (comma ? ",\n" : "") << "    {\"v\":" << v << ",\"xy\":[" << v%n << "," << v/n << "]";
    std::cout << ",\"B_v_from_D\":\"" << dec(r.B) << "\",\"A2_v_from_fibres\":\"" << dec(r.A) << "\",\"histogram_L\":[";
    for (int t=1;t<(int)r.hist.size();++t) { if(t>1) std::cout << ','; std::cout << '"' << dec(r.hist[t]) << '"'; }
    std::cout << "],\"root_B\":[";
    bool first=true; for(int p=0;p<(int)r.rootB.size();++p) if(r.rootB[p]) {
        if(!first) std::cout << ','; first=false; std::cout << "[" << p << ",\"" << dec(r.rootB[p]) << "\"]";
    }
    std::cout << "]}";
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: tail_full_dp n v|all\n"; return 2; }
    try {
        int n=std::stoi(argv[1]); std::vector<int> holes;
        if(std::string(argv[2])=="all") for(int v=0;v<n*n;++v) holes.push_back(v);
        else holes.push_back(std::stoi(argv[2]));
        std::cout << "{\"n\":" << n << ",\"method\":\"shared-f suffix DP plus compressed pre-neighbour prefix DP\",\"u128_checked\":true,\"holes\":[";
        bool comma=false; for(int v:holes) { auto r=compute_hole(n,v); emit(r,n,comma); comma=true; }
        std::cout << "]}\n";
    } catch(const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 1; }
}
