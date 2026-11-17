// Exact backward subset DP for the two-row-release completion formula.
// Vertices in H are indexed by (x,y), 0<=x<n, 2<=y<n; M=n(n-2).
// W[p-1][S] stores the weighted completion sum from an already-built S
// having p births (p=1,2,3).  U[S] stores the unweighted completion count.
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using u128 = unsigned __int128;
constexpr int PMAX = 3;

static std::string to_string(u128 x) {
  if (!x) return "0";
  std::string s;
  while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
  std::reverse(s.begin(), s.end()); return s;
}
static u128 fact(int n) {
  u128 r=1; for (int i=2;i<=n;i++) r*=u128(i); return r;
}
static const u128 LIMIT30 = fact(30); // All requested counts are <30!.
static u128 add_checked(u128 a, u128 b) {
  if (b > LIMIT30 || a > std::numeric_limits<u128>::max()-b || a > LIMIT30-b)
    throw std::overflow_error("u128 addition exceeds 30! bound");
  return a+b;
}
static u128 mul_checked(u128 a, uint32_t b) {
  if (b && (a > std::numeric_limits<u128>::max()/b || a > LIMIT30/b))
    throw std::overflow_error("u128 multiplication exceeds 30! bound");
  return a*b;
}

struct Problem {
  int n, M, maxp;
  std::vector<uint32_t> hnb;
  std::vector<uint8_t> rnb;
  explicit Problem(int q):n(q),M(q*(q-2)),maxp((q+1)/2),hnb(M),rnb(M) {
    if (n<3 || n>6 || M>24 || maxp>PMAX) throw std::invalid_argument("n outside supported range");
    auto ix=[&](int x,int y){ return (y-2)*n+x; };
    for(int y=2;y<n;y++) for(int x=0;x<n;x++) {
      int i=ix(x,y); rnb[i]=(y==2);
      for(auto d: std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}) {
        int X=x+d[0], Y=y+d[1];
        if (0<=X && X<n && 2<=Y && Y<n) hnb[i] |= uint32_t(1)<<ix(X,Y);
      }
    }
  }
};

static std::vector<uint32_t> masks_of_size(int M,int k) {
  std::vector<uint32_t> out;
  if(k==0) { out.push_back(0); return out; }
  uint32_t comb=(uint32_t(1)<<k)-1, lim=uint32_t(1)<<M;
  while(comb<lim) {
    out.push_back(comb);
    uint32_t c=comb&-comb, r=comb+c;
    if(r>=lim) break;
    comb=(((r^comb)>>2)/c)|r;
  }
  return out;
}

struct Result { u128 HR=0, HT=0; size_t states=0; };
static Result run(const Problem& g, int threads) {
  const size_t Z=size_t(1)<<g.M; const uint32_t full=uint32_t(Z-1);
  // n=6 would allocate 3*256MiB + 256MiB plus small layer buffers.
  std::array<std::vector<u128>,PMAX> W;
  for(auto& a:W) a.assign(Z,0);
  std::vector<u128> U(Z,0);
  for(int p=1;p<=g.maxp;p++) W[p-1][full]=fact(g.M+g.n)/fact(g.M+p);
  U[full]=1;
  std::atomic<bool> bad{false}; std::exception_ptr ep; std::mutex ep_mu;
  size_t states=1;
  for(int k=g.M-1;k>=1;k--) {
    auto layer=masks_of_size(g.M,k); states+=layer.size();
    auto worker=[&](size_t lo,size_t hi) {
      try {
        for(size_t z=lo;z<hi && !bad.load();z++) {
          uint32_t S=layer[z]; u128 us=0; std::array<u128,PMAX> ws{};
          uint32_t absent=full^S;
          while(absent) {
            uint32_t bit=absent&-absent; absent^=bit;
            int i=__builtin_ctz(bit); bool allowed=g.rnb[i] || (g.hnb[i]&S);
            if(!allowed) continue;
            uint32_t T=S|bit; us=add_checked(us,U[T]);
            bool birth=(g.hnb[i]&S)==0;
            for(int p=1;p<=g.maxp;p++) {
              if(birth) {
                if(p<g.maxp) ws[p-1]=add_checked(ws[p-1],mul_checked(W[p][T],uint32_t(k+p+1)));
              } else ws[p-1]=add_checked(ws[p-1],W[p-1][T]);
            }
          }
          U[S]=us; for(int p=1;p<=g.maxp;p++) W[p-1][S]=ws[p-1];
        }
      } catch(...) { bad.store(true); std::lock_guard<std::mutex> lk(ep_mu); if(!ep)ep=std::current_exception(); }
    };
    int t=std::min<int>(threads,layer.size()); std::vector<std::thread> pool;
    for(int j=0;j<t;j++) { size_t lo=layer.size()*j/t, hi=layer.size()*(j+1)/t; pool.emplace_back(worker,lo,hi); }
    for(auto& th:pool) th.join(); if(ep) std::rethrow_exception(ep);
  }
  // Empty state has p=0. Its only permitted first additions are births from R.
  Result r; r.states=states+1;
  for(int i=0;i<g.M;i++) if(g.rnb[i]) { uint32_t bit=uint32_t(1)<<i; r.HR=add_checked(r.HR,U[bit]); r.HT=add_checked(r.HT,W[0][bit]); }
  return r;
}

int main(int argc,char**argv) {
  if(argc!=5) { std::cerr<<"usage: row_release_native n threads output.json run|dry-run\n"; return 2; }
  int n=std::stoi(argv[1]), T=std::stoi(argv[2]); std::string out=argv[3], mode=argv[4];
  if(T<1||T>64) return 2; Problem g(n);
  if(mode=="dry-run") {
    size_t states=size_t(1)<<g.M, bytes=(size_t(g.maxp)+1)*states*sizeof(u128);
    std::ofstream f(out); f<<"{\"n\":"<<n<<",\"M\":"<<g.M<<",\"maxp\":"<<g.maxp<<",\"array_bytes\":"<<bytes<<",\"limit30\":\""<<to_string(LIMIT30)<<"\"}\n"; return 0;
  }
  if(mode!="run") return 2;
  Result r=run(g,T); std::ofstream f(out);
  f<<"{\"n\":"<<n<<",\"M\":"<<g.M<<",\"threads\":"<<T<<",\"algorithm\":\"backward-subset-u128\",\"states_visited\":"<<r.states<<",\"H_two_rows\":\""<<to_string(r.HR)<<"\",\"H_toprow\":\""<<to_string(r.HT)<<"\",\"limit30\":\""<<to_string(LIMIT30)<<"\"}\n";
}
