// Exact peak-set inclusion-exclusion through independent sets.
// For n<=6 the largest table is the 5,598,861 independent sets of 6x6.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct U192 { std::array<uint64_t,3> x{}; };
struct U256 { std::array<uint64_t,4> x{}; };

static U192 add192(U192 a, const U192& b) {
  for (int i=0;i<3;i++) { uint64_t old=a.x[i]; a.x[i]+=b.x[i]; if (a.x[i]<old) {
      for (int j=i+1;j<3;j++) { ++a.x[j]; if (a.x[j]) break; }
    }} return a;
}
static U192 mulsmall192(U192 a, uint64_t q) {
  U192 r; __uint128_t c=0; for(int i=0;i<3;i++) { __uint128_t z=(__uint128_t)a.x[i]*q+c; r.x[i]=(uint64_t)z;c=z>>64; }
  assert(c==0); return r;
}
static U192 divsmall192(U192 a, uint64_t q) {
  U192 r; __uint128_t rem=0; for(int i=2;i>=0;i--) { __uint128_t z=(rem<<64)|a.x[i]; r.x[i]=(uint64_t)(z/q); rem=z%q; }
  assert(rem==0); return r;
}
static bool iszero(const U192&a){return a.x[0]==0&&a.x[1]==0&&a.x[2]==0;}
static void addscaled(U256& a,const U192&b,uint64_t q) {
  __uint128_t carry=0; for(int i=0;i<3;i++){ __uint128_t z=(__uint128_t)b.x[i]*q+carry; uint64_t lo=(uint64_t)z; carry=z>>64; uint64_t old=a.x[i];a.x[i]+=lo;if(a.x[i]<old)++carry; }
  uint64_t old=a.x[3]; a.x[3]+=(uint64_t)carry; assert(a.x[3]>=old);
}
static int cmp(const U256&a,const U256&b){for(int i=3;i>=0;i--)if(a.x[i]!=b.x[i])return a.x[i]<b.x[i]?-1:1;return 0;}
static U256 sub(U256 a,const U256&b){uint64_t borrow=0;for(int i=0;i<4;i++){uint64_t bi=b.x[i]+borrow;bool o=bi<b.x[i];uint64_t old=a.x[i];a.x[i]-=bi;borrow=o||old<bi;}assert(!borrow);return a;}
static U256 add256(U256 a,const U256&b){for(int i=0;i<4;i++){uint64_t old=a.x[i];a.x[i]+=b.x[i];if(a.x[i]<old){for(int j=i+1;j<4;j++){++a.x[j];if(a.x[j])break;}}}return a;}
static U256 divsmall256(U256 a,uint64_t q) { U256 r; __uint128_t rem=0; for(int i=3;i>=0;i--) { __uint128_t z=(rem<<64)|a.x[i]; r.x[i]=(uint64_t)(z/q); rem=z%q; } assert(rem==0); return r; }
static std::string str(U256 a){if(a.x==std::array<uint64_t,4>{})return"0";std::string s;while(a.x!=std::array<uint64_t,4>{}){__uint128_t r=0;for(int i=3;i>=0;i--){__uint128_t z=(r<<64)|a.x[i];a.x[i]=(uint64_t)(z/10);r=z%10;}s.push_back(char('0'+(int)r));}std::reverse(s.begin(),s.end());return s;}
static std::string signedstr(const U256&p,const U256&m){int c=cmp(p,m);return c<0?"-"+str(sub(m,p)):str(sub(p,m));}
static U256 widen(const U192&a){U256 b;b.x[0]=a.x[0];b.x[1]=a.x[1];b.x[2]=a.x[2];return b;}

static int n,N,alpha;
static std::vector<uint64_t> valid;
static std::vector<std::vector<uint64_t>> bags;
static void enumerate(int row,uint64_t previous,uint64_t set){
  if(row==n){bags[__builtin_popcountll(set)].push_back(set);return;}
  for(uint64_t r:valid)if(!(r&previous))enumerate(row+1,r,set|(r<<(row*n)));
}
static uint64_t neighbors(uint64_t s){
  uint64_t col0=0,colLast=0;for(int y=0;y<n;y++){col0|=1ULL<<(y*n);colLast|=1ULL<<(y*n+n-1);}
  uint64_t ans=((s&~col0)>>1)|((s&~colLast)<<1);
  ans|=(s<<n);ans|=(s>>n);return ans&((1ULL<<N)-1);
}
static uint64_t choose(int a,int b){if(b<0||b>a)return 0;b=std::min(b,a-b);uint64_t z=1;for(int i=1;i<=b;i++)z=z*(a-b+i)/i;return z;}
int main(int argc,char**argv){
  if(argc!=2)return 2;n=std::atoi(argv[1]);if(n<2||n>6)return 2;N=n*n;alpha=(N+1)/2;
  for(uint64_t r=0;r<(1ULL<<n);r++)if(!(r&(r<<1)))valid.push_back(r);
  bags.resize(alpha+1);enumerate(0,0,0);
  size_t total=0;for(auto&b:bags)total+=b.size();
  std::vector<uint64_t> sets;sets.reserve(total);std::vector<size_t> start(alpha+2);for(int k=0;k<=alpha;k++){start[k]=sets.size();sets.insert(sets.end(),bags[k].begin(),bags[k].end());std::vector<uint64_t>().swap(bags[k]);}start[alpha+1]=sets.size();
  std::unordered_map<uint64_t,uint32_t> index;index.max_load_factor(.7);index.reserve(total*2);for(uint32_t i=0;i<sets.size();i++)index.emplace(sets[i],i);
  std::vector<U192> C(total);C[0].x[0]=1;for(int i=2;i<=N;i++)C[0]=mulsmall192(C[0],i);
  for(int k=1;k<=alpha;k++)for(size_t z=start[k];z<start[k+1];z++){
    uint64_t s=sets[z],q=s;U192 sum;while(q){uint64_t bit=q&-q;q-=bit;sum=add192(sum,C[index.at(s^bit)]);}int den=k+__builtin_popcountll(neighbors(s));C[z]=divsmall192(sum,den);
  }
  // Per-puncture inclusion-exclusion.  A: exactly one secondary peak after
  // the prescribed top label v. B: the prescribed rank-2 label is outside
  // N[v], and there are no further peaks.  Positive/negative streams avoid
  // signed-integer overflow before the exact final cancellation.
  uint64_t full=(1ULL<<N)-1;
  std::vector<uint64_t> closedv(N);
  for(int v=0;v<N;v++) closedv[v]=(1ULL<<v)|neighbors(1ULL<<v);
  std::vector<U256> ap(N),am(N),bp(N),bm(N);
  for(size_t z=0;z<sets.size();z++) {
    uint64_t s=sets[z], cl=s|neighbors(s); int k=__builtin_popcountll(s);
    bool aneg=(k%2==0);          // (-1)^(k-1); the k=0 term is zero.
    bool bneg=(k%2==1);          // (-1)^k.
    for(int v=0;v<N;v++) if(!(s&closedv[v])) {
      if(k) addscaled(aneg?am[v]:ap[v],C[z],k);
      uint64_t avail=full&~(closedv[v]|cl);
      addscaled(bneg?bm[v]:bp[v],C[z],__builtin_popcountll(avail));
    }
  }
  std::vector<U256> av(N),bv(N); U256 atot,btot;
  for(int v=0;v<N;v++) {
    assert(cmp(ap[v],am[v])>=0 && cmp(bp[v],bm[v])>=0);
    av[v]=divsmall256(sub(ap[v],am[v]),N);
    bv[v]=divsmall256(sub(bp[v],bm[v]),uint64_t(N)*(N-1));
    atot=add256(atot,av[v]); btot=add256(btot,bv[v]);
  }
  std::cout<<"{\"n\":"<<n<<",\"N\":"<<N<<",\"independent_sets\":"<<total<<",\"A2\":\""<<str(atot)<<"\",\"B\":\""<<str(btot)<<"\",\"per_vertex\":[";
  for(int v=0;v<N;v++) { if(v) std::cout<<','; std::cout<<"{\"v\":"<<v<<",\"A2_v\":\""<<str(av[v])<<"\",\"B_v\":\""<<str(bv[v])<<"\"}"; }
  std::cout<<"]}\n";
  return 0;
  std::vector<U256> plus(alpha+1),minus(alpha+1),rp(N),rm(N);int diam=2*(n-1);std::vector<U256> dp(diam+1),dm(diam+1);U256 sqp,sqm;
  for(size_t z=0;z<sets.size();z++){
    int k=__builtin_popcountll(sets[z]); if(iszero(C[z]))continue;
    for(int r=1;r<=k;r++){auto&dst=((k-r)&1)?minus[r]:plus[r];addscaled(dst,C[z],choose(k,r));}
    std::array<int,18> vs{};int nv=0;uint64_t q=sets[z];while(q){int v=__builtin_ctzll(q);q&=q-1;vs[nv++]=v;}
    bool rootneg=((k-1)&1); // (-1)^(k-1)
    for(int i=0;i<nv;i++)addscaled(rootneg?rm[vs[i]]:rp[vs[i]],C[z],1);
    if(k<2)continue;
    bool neg=(k&1); // (-1)^(k-2)
    uint64_t squares=0;
    for(int i=0;i<nv;i++)for(int j=i+1;j<nv;j++){
      int d=std::abs(vs[i]%n-vs[j]%n)+std::abs(vs[i]/n-vs[j]/n);
      addscaled(neg?dm[d]:dp[d],C[z],1);
      int dx=vs[i]%n-vs[j]%n,dy=vs[i]/n-vs[j]/n;squares+=uint64_t(dx*dx+dy*dy);
    }
    addscaled(neg?sqm:sqp,C[z],squares);
  }
  uint64_t corner=1, two=corner|(1ULL<<2), nine=0, even=0, odd=0;
  for(int y=0;y<n;y++)for(int x=0;x<n;x++){uint64_t bit=1ULL<<(y*n+x);if((x+y)&1)odd|=bit;else even|=bit;if(n==6&&x%2==0&&y%2==0)nine|=bit;}
  std::array<std::pair<const char*,uint64_t>,5> probes{{{"corner",corner},{"two_same_row",two},{"even_even_3x3",nine},{"checkerboard_even",even},{"checkerboard_odd",odd}}};
  std::cout<<"{\"n\":"<<n<<",\"N\":"<<N<<",\"independent_sets\":"<<total<<",\"peak_histogram\":{";
  for(int k=1;k<=alpha;k++){if(k>1)std::cout<<',';std::cout<<'\"'<<k<<"\":\""<<signedstr(plus[k],minus[k])<<"\"";}
  std::cout<<"},\"A1\":\""<<signedstr(plus[1],minus[1])<<"\",\"A2\":\""<<signedstr(plus[2],minus[2])<<"\",\"distance_counts\":{";
  for(int d=0;d<=diam;d++)if(cmp(dp[d],dm[d])||d==diam){static bool first=true;if(!first)std::cout<<',';first=false;std::cout<<'\"'<<d<<"\":\""<<signedstr(dp[d],dm[d])<<"\"";}
  std::cout<<"},\"squared_euclidean_total\":\""<<signedstr(sqp,sqm)<<"\",\"single_peak_root_counts\":[";
  for(int v=0;v<N;v++){if(v)std::cout<<',';std::cout<<'\"'<<signedstr(rp[v],rm[v])<<'\"';}
  std::cout<<"],\"selected_joint_peak_counts\":{";
  bool first=true;for(auto [name,s]:probes)if(s&&__builtin_popcountll(s)<=alpha&&!(s&neighbors(s))){if(!first)std::cout<<',';first=false;std::cout<<'\"'<<name<<"\":\""<<str(widen(C[index.at(s)]))<<"\"";}
  std::cout<<"}}\n";
}
