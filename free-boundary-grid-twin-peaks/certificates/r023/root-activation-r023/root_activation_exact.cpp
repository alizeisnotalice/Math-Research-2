// Exact h-completion scan for isolated low roots on P_n square, n=3..5.
#include <array>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using U128 = unsigned __int128;

static std::string dec(U128 x) {
  if (!x) return "0";
  std::string s;
  while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
  return {s.rbegin(), s.rend()};
}

// Exact non-negative products up to 512 bits, used only for fraction ordering.
// Base 2^32 makes every elementary product-plus-carry fit in U128.
struct Wide {
  std::array<uint32_t,16> w{};
  explicit Wide(U128 x=0) { for(int i=0;i<4;i++){w[i]=uint32_t(x);x>>=32;} }
  void mul(U128 x) {
    std::array<uint32_t,4> b{};for(int i=0;i<4;i++){b[i]=uint32_t(x);x>>=32;}
    std::array<uint32_t,16> r{};
    for(int i=0;i<16;i++) if(w[i]) {
      U128 carry=0;
      for(int j=0;j<4&&i+j<16;j++) {
        U128 z=U128(r[i+j])+U128(w[i])*b[j]+carry;
        r[i+j]=uint32_t(z); carry=z>>32;
      }
      int at=i+4;
      while(carry&&at<16) { U128 z=U128(r[at])+carry; r[at]=uint32_t(z);carry=z>>32;at++; }
      if(carry){std::cerr<<"wide overflow\n";std::abort();}
    }
    w=r;
  }
  friend bool operator<(const Wide&a,const Wide&b) {
    for(int i=15;i>=0;i--)if(a.w[i]!=b.w[i])return a.w[i]<b.w[i];return false;
  }
};
static Wide prod(std::initializer_list<U128> x){Wide z(1);for(U128 a:x)z.mul(a);return z;}

struct Rec { bool set=false; uint32_t C=0; int y=0; U128 num=0, den=1; };

struct Grid {
  int n, N; uint32_t full, col0=0, coll=0; std::vector<uint32_t> adj;
  explicit Grid(int q):n(q),N(q*q),full((uint32_t(1)<<N)-1),adj(N) {
    for (int yy=0;yy<n;yy++) for (int xx=0;xx<n;xx++) {
      int v=yy*n+xx;
      for (auto d: {std::pair{1,0},std::pair{-1,0},std::pair{0,1},std::pair{0,-1}}) {
        int X=xx+d.first,Y=yy+d.second;
        if (X>=0&&X<n&&Y>=0&&Y<n) adj[v]|=uint32_t(1)<<(Y*n+X);
      }
      if (xx==0) col0|=uint32_t(1)<<v;
      if (xx==n-1) coll|=uint32_t(1)<<v;
    }
  }
  uint32_t nbr(uint32_t s) const {
    return (((s&~col0)>>1)|((s&~coll)<<1)|(s<<n)|(s>>n))&full;
  }
  bool connected(uint32_t s) const {
    uint32_t x=s&-s;
    while (true) {
      uint32_t y=x|(nbr(x)&s);
      if (y==x) return x==s;
      x=y;
    }
  }
};

static bool less_p(const Rec&a,const Rec&b) {
  return prod({a.num,b.den}) < prod({b.num,a.den});
}
static bool less_linear(const Rec&a,const Rec&b,int ka,int kb) {
  return prod({U128(ka+1),a.num,b.den}) < prod({U128(kb+1),b.num,a.den});
}
static bool less_sqrt(const Rec&a,const Rec&b,int ka,int kb) {
  // Compare k_a p_a^2 and k_b p_b^2 exactly, without square roots.
  return prod({U128(ka),a.num,a.num,b.den,b.den}) < prod({U128(kb),b.num,b.num,a.den,a.den});
}

static void update(Rec& dst,const Rec& cand,int k,int which) {
  if (!dst.set || (which==0?less_p(cand,dst):which==1?less_linear(cand,dst,k,k):less_sqrt(cand,dst,k,k))) dst=cand;
}

static std::string frac(U128 x,U128 y) { return dec(x)+"/"+dec(y); }
static void coords(std::ostream&o,uint32_t C,const Grid&g) {
  o<<'['; bool first=true; for(int v=0;v<g.N;v++) if(C&(uint32_t(1)<<v)) {
    if(!first)o<<','; first=false; o<<'['<<(v%g.n)<<','<<(v/g.n)<<']';
  } o<<']';
}

static void emit_rec(std::ostream&o,const Rec&r,int k,const Grid&g,const std::vector<U128>&h) {
  uint32_t ybit=uint32_t(1)<<r.y,S=r.C|ybit;
  o<<"{\"C_coordinates\":"; coords(o,r.C,g);
  o<<",\"y\":["<<(r.y%g.n)<<','<<(r.y/g.n)<<"]"
   <<",\"C_size\":"<<k<<",\"h_S\":\""<<dec(r.den)<<"\""
   <<",\"activation_numerator\":\""<<dec(r.num)<<"\""
   <<",\"p_act\":\""<<frac(r.num,r.den)<<"\""
   <<",\"p_times_C_plus_1\":\""<<frac(U128(k+1)*r.num,r.den)<<"\""
   <<",\"p_squared_comparison_exact\":true"
   <<",\"neighbour_completion_terms\":[";
  bool first=true; uint32_t zset=g.adj[r.y];
  while(zset) { uint32_t bit=zset&-zset; zset^=bit; int z=__builtin_ctz(bit);
    if(!first)o<<','; first=false;
    o<<"{\"z\":["<<(z%g.n)<<','<<(z/g.n)<<"],\"h_S_plus_z\":\""<<dec(h[S|bit])
     <<"\",\"also_adjacent_C\":"<<((g.adj[z]&r.C)?"true":"false")<<'}';
  }
  o<<"]}";
}

int main(int argc,char**argv) {
  if(argc!=3) return 2; int n=std::atoi(argv[1]); if(n<3||n>5)return 2;
  Grid g(n); const uint64_t cap=uint64_t(g.full)+1;
  std::vector<U128> h(cap); h[g.full]=1;
  for(uint32_t s=g.full-1;s;s--) {
    uint32_t b=g.nbr(s)&~s; U128 z=0;
    while(b){uint32_t bit=b&-b;b^=bit;z+=h[s|bit];}
    h[s]=z;
  }
  std::vector<Rec> mp(g.N+1), ml(g.N+1), ms(g.N+1);
  std::vector<uint64_t> cnt(g.N+1); uint64_t checked=0;
  for(int y=0;y<g.N;y++) {
    uint32_t avail=g.full ^ ((uint32_t(1)<<y)|g.adj[y]);
    for(uint32_t C=avail;C;C=(C-1)&avail) {
      if(!g.connected(C)) continue;
      int k=__builtin_popcount(C); uint32_t S=C|(uint32_t(1)<<y);
      U128 num=0; uint32_t zset=g.adj[y];
      while(zset){uint32_t bit=zset&-zset;zset^=bit;num+=h[S|bit];}
      Rec r{true,C,y,num,h[S]};
      if(!(r.num>0&&r.num<=r.den)){std::cerr<<"bad probability\n";return 3;}
      update(mp[k],r,k,0); update(ml[k],r,k,1); update(ms[k],r,k,2); cnt[k]++; checked++;
    }
  }
  std::ofstream o(argv[2]);
  o<<"{\"n\":"<<n<<",\"N\":"<<g.N
   <<",\"method\":\"exact U128 h table plus isolated-root scan; p comparisons use U256/U512\""
   <<",\"h_table_bytes\":"<<cap*sizeof(U128)<<",\"eligible_pairs\":"<<checked<<",\"by_C_size\":[";
  bool first=true; for(int k=1;k<=g.N;k++) if(cnt[k]){
    if(!first)o<<',';first=false; o<<"{\"C_size\":"<<k<<",\"eligible_pairs\":"<<cnt[k]
      <<",\"minimum_p_act\":";emit_rec(o,mp[k],k,g,h);
    o<<",\"minimum_p_times_C_plus_1\":";emit_rec(o,ml[k],k,g,h);
    o<<",\"minimum_p_times_sqrt_C\":";emit_rec(o,ms[k],k,g,h);o<<'}';
  }
  o<<"],\"checks\":{\"connected_C_only\":true,\"isolated_y_only\":true,\"all_probabilities_in_0_1\":true}}\n";
}
