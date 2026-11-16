// Exact vertical-cut covariance, reduced by the D4 action.
// Every stored mask is the least member of its D4 orbit.  C(I)=49! P(I all
// peaks) is computed cardinality layer by layer from C(I\{v}).  For a
// non-D4-invariant cut statistic, every one of the eight transforms is
// accumulated; a representative is never merely multiplied by its orbit size.
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using Mask = uint64_t;
constexpr int MAXN=49, MAXALPHA=25, THREADS=8;
static int n=7, N=49, ALPHA=25;
constexpr Mask EMPTY = ~Mask(0);

struct U256 { std::array<uint64_t,4> w{}; };
static bool zero(const U256&a){return !(a.w[0]|a.w[1]|a.w[2]|a.w[3]);}
static void add(U256&a,const U256&b){
  uint64_t carry=0;
  for(int i=0;i<4;i++){ uint64_t x=a.w[i]; uint64_t y=b.w[i]; uint64_t z=x+y; bool c1=z<x; uint64_t q=z+carry; bool c2=q<z; a.w[i]=q; carry=c1||c2; }
  assert(!carry);
}
static U256 mulsmall(const U256&a,uint64_t q){
  U256 r; __uint128_t c=0;
  for(int i=0;i<4;i++){ __uint128_t z=(__uint128_t)a.w[i]*q+c; r.w[i]=(uint64_t)z;c=z>>64; }
  assert(!c); return r;
}
static U256 divsmall(const U256&a,uint64_t q){
  U256 r; __uint128_t rem=0;
  for(int i=3;i>=0;i--){ __uint128_t z=(rem<<64)|a.w[i];r.w[i]=(uint64_t)(z/q);rem=z%q; }
  assert(!rem);return r;
}
static std::string str(U256 a){
  if(zero(a))return "0";std::string s;
  while(!zero(a)){__uint128_t r=0;for(int i=3;i>=0;i--){__uint128_t z=(r<<64)|a.w[i];a.w[i]=(uint64_t)(z/10);r=z%10;}s.push_back(char('0'+(int)r));}
  std::reverse(s.begin(),s.end());return s;
}
static int cmp(const U256&a,const U256&b){for(int i=3;i>=0;i--)if(a.w[i]!=b.w[i])return a.w[i]<b.w[i]?-1:1;return 0;}
static U256 sub(U256 a,const U256&b){uint64_t borrow=0;for(int i=0;i<4;i++){uint64_t bi=b.w[i]+borrow;bool c=bi<b.w[i];uint64_t old=a.w[i];a.w[i]-=bi;borrow=c||(old<bi);}assert(!borrow);return a;}
static std::string sstr(const U256&p,const U256&m){return cmp(p,m)<0?"-"+str(sub(m,p)):str(sub(p,m));}
static void addsigned(U256& p,U256& m,const U256& v,int64_t q){
  if(q>=0) add(p,mulsmall(v,uint64_t(q))); else add(m,mulsmall(v,uint64_t(-q)));
}

struct U512 { std::array<uint64_t,8> w{}; };
static void addword(U512& a,int at,uint64_t v){while(v){assert(at<8);uint64_t old=a.w[at];a.w[at]+=v;v=(a.w[at]<old);++at;}}
static U512 mul256(const U256&a,const U256&b){U512 z;for(int i=0;i<4;i++)for(int j=0;j<4;j++){__uint128_t q=(__uint128_t)a.w[i]*b.w[j];addword(z,i+j,uint64_t(q));addword(z,i+j+1,uint64_t(q>>64));}return z;}
static U512 twice(U512 a){U512 z;for(int i=0;i<8;i++){addword(z,i,a.w[i]);addword(z,i,a.w[i]);}return z;}
static int cmp512(const U512&a,const U512&b){for(int i=7;i>=0;i--)if(a.w[i]!=b.w[i])return a.w[i]<b.w[i]?-1:1;return 0;}
static U512 sub512(U512 a,const U512&b){uint64_t borrow=0;for(int i=0;i<8;i++){uint64_t bi=b.w[i]+borrow;bool over=bi<b.w[i];uint64_t old=a.w[i];a.w[i]-=bi;borrow=over||old<bi;}assert(!borrow);return a;}
static std::string str512(U512 a){if(a.w==std::array<uint64_t,8>{})return "0";std::string s;while(a.w!=std::array<uint64_t,8>{}){__uint128_t r=0;for(int i=7;i>=0;i--){__uint128_t q=(r<<64)|a.w[i];a.w[i]=uint64_t(q/10);r=q%10;}s.push_back(char('0'+int(r)));}std::reverse(s.begin(),s.end());return s;}
struct SM256 { int sign; U256 mag; };
struct SM512 { int sign; U512 mag; };
static SM256 signmag(const U256&p,const U256&m){int c=cmp(p,m);return {c,c<0?sub(m,p):sub(p,m)};}
static SM512 product(SM256 a,SM256 b){return {a.sign*b.sign,mul256(a.mag,b.mag)};}
static SM512 subtract_signed(SM512 a,SM512 b){b.sign=-b.sign;if(!a.sign)return b;if(!b.sign)return a;if(a.sign==b.sign){U512 z=a.mag;for(int i=0;i<8;i++)addword(z,i,b.mag.w[i]);return {a.sign,z};}int c=cmp512(a.mag,b.mag);if(!c)return {0,{}};return c>0?SM512{a.sign,sub512(a.mag,b.mag)}:SM512{b.sign,sub512(b.mag,a.mag)};}
static std::string sstr512(SM512 z){return z.sign<0?"-"+str512(z.mag):str512(z.mag);}

static std::array<Mask,MAXN> closedv;
static std::array<Mask,MAXN> mapbit[8];
static std::vector<Mask> rows;
static uint64_t raw_by_k[MAXALPHA+1]={
  1,49,1092,14690,133544,870589,4216498,15516804,44031035,
  97284860,168434824,229465420,246654832,209621800,141138510,
  75497964,32229088,11059218,3084524,710124,137368,22611,3134,344,26,1};

static Mask transform(Mask s,int op){Mask r=0;while(s){Mask b=s&-s;s-=b;int v=__builtin_ctzll(b);r|=mapbit[op][v];}return r;}
static Mask canon(Mask s,uint8_t*orb=nullptr){
  Mask vals[8];vals[0]=s;Mask mn=s;
  for(int op=1;op<8;op++){vals[op]=transform(s,op);if(vals[op]<mn)mn=vals[op];}
  if(orb){*orb=0;for(int i=0;i<8;i++){bool fresh=true;for(int j=0;j<i;j++)if(vals[j]==vals[i]){fresh=false;break;}if(fresh)(*orb)++;}}
  return mn;
}

struct Layer { std::vector<Mask> m; std::vector<uint8_t> orb; };
static std::array<Layer,MAXALPHA+1> layers;
static std::atomic<int> next_task{0};
static std::vector<std::pair<Mask,Mask>> starts;
static uint64_t hash64(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;return x;}

static void enum_tail(int row, Mask previous, Mask set, std::array<std::vector<Mask>,MAXALPHA+1>&localm,
                      std::array<std::vector<uint8_t>,MAXALPHA+1>&localo){
  if(row==n){uint8_t o;if(canon(set,&o)==set){int k=__builtin_popcountll(set);localm[k].push_back(set);localo[k].push_back(o);}return;}
  for(Mask r:rows)if(!(r&previous))enum_tail(row+1,r,set|(r<<(row*n)),localm,localo);
}
static void enumerate_worker(){
  std::array<std::vector<Mask>,MAXALPHA+1> lm;
  std::array<std::vector<uint8_t>,MAXALPHA+1> lo;
  for(;;){int t=next_task.fetch_add(1);if(t>=(int)starts.size())break;auto [a,b]=starts[t];enum_tail(2,b,a|(b<<n),lm,lo);}
  static std::mutex mu;
  std::lock_guard<std::mutex> g(mu);
  for(int k=0;k<=ALPHA;k++){
    layers[k].m.insert(layers[k].m.end(),lm[k].begin(),lm[k].end());
    layers[k].orb.insert(layers[k].orb.end(),lo[k].begin(),lo[k].end());
  }
}

// The immutable previous-layer masks carry hash keys; slots store index+1.
struct Table{
  std::vector<std::atomic<uint32_t>> slot;const std::vector<Mask>* keys=nullptr;uint64_t mask=0;
  void build(const std::vector<Mask>&a){
    keys=&a;uint64_t cap=1;while(cap<a.size()*2)cap<<=1;slot=std::vector<std::atomic<uint32_t>>(cap);mask=cap-1;
    auto parallel=[&](size_t total,auto body){std::atomic<size_t> next{0};std::vector<std::thread> ts;for(int q=0;q<THREADS;q++)ts.emplace_back([&]{for(;;){size_t lo=next.fetch_add(1<<16);if(lo>=total)break;body(lo,std::min(total,lo+(1<<16)));}});for(auto&t:ts)t.join();};
    parallel(cap,[&](size_t lo,size_t hi){for(size_t i=lo;i<hi;i++)slot[i].store(0,std::memory_order_relaxed);});
    parallel(a.size(),[&](size_t lo,size_t hi){for(size_t i=lo;i<hi;i++){uint64_t h=hash64(a[i])&mask;uint32_t wanted=uint32_t(i+1);for(;;){uint32_t old=slot[h].load(std::memory_order_relaxed);if(!old){if(slot[h].compare_exchange_weak(old,wanted,std::memory_order_relaxed,std::memory_order_relaxed))break;}else{assert((*keys)[old-1]!=a[i]);h=(h+1)&mask;}}}});
  }
  uint32_t get(Mask x)const{uint64_t h=hash64(x)&mask;for(;;){uint32_t p=slot[h].load(std::memory_order_relaxed);assert(p);if((*keys)[p-1]==x)return p-1;h=(h+1)&mask;}}
};

static Mask neighbourhood(Mask s){Mask z=0;while(s){Mask b=s&-s;s-=b;z|=closedv[__builtin_ctzll(b)];}return z;}
static uint64_t choose(int a,int b){if(b<0||b>a)return 0;b=std::min(b,a-b);uint64_t x=1;for(int i=1;i<=b;i++)x=x*(a-b+i)/i;return x;}
constexpr int MAXCUT=6;
struct Accum { U256 a2p{},a2m{};std::array<U256,MAXCUT> lp{},lm{};std::array<std::array<U256,MAXCUT>,MAXCUT> hp{},hm{}; };

int main(int argc,char**argv){
  if(argc!=2)return 2; n=std::atoi(argv[1]); if(n<3||n>7)return 2;
  N=n*n; ALPHA=(N+1)/2;
  for(int y=0;y<n;y++)for(int x=0;x<n;x++){
    int v=y*n+x;closedv[v]=Mask(1)<<v;
    for(auto [X,Y]:{std::pair{x-1,y},{x+1,y},{x,y-1},{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)closedv[v]|=Mask(1)<<(Y*n+X);
    for(int op=0;op<8;op++){int X=x,Y=y;if(op>=4)X=n-1-X;for(int q=0;q<op%4;q++){int T=n-1-Y;Y=X;X=T;}mapbit[op][v]=Mask(1)<<(Y*n+X);}
  }
  for(Mask r=0;r<(1<<n);r++)if(!(r&(r<<1)))rows.push_back(r);
  for(Mask a:rows)for(Mask b:rows)if(!(a&b))starts.push_back({a,b});
  std::cerr<<"enumerating "<<starts.size()<<" two-row tasks with "<<THREADS<<" threads\n";
  std::vector<std::thread> th;for(int i=0;i<THREADS;i++)th.emplace_back(enumerate_worker);for(auto&t:th)t.join();
  uint64_t weighted=0;
  for(int k=0;k<=ALPHA;k++){uint64_t got=0;for(auto o:layers[k].orb)got+=o;weighted+=got;if(n==7)assert(got==raw_by_k[k]);std::cerr<<"k="<<k<<" reps="<<layers[k].m.size()<<" weighted="<<got<<"\n";}
  if(n==7)assert(weighted==1280128950ULL);
  std::vector<U256> prev(1),cur;
  U256 fact;fact.w[0]=1;for(int i=2;i<=N;i++)fact=mulsmall(fact,i);prev[0]=fact;
  U256 a2p{},a2m{};
  std::array<U256,MAXCUT> lp{},lm{};
  std::array<std::array<U256,MAXCUT>,MAXCUT> hp{},hm{};
  // k=0 contributes neither peaks nor pairs.  All other layers are independent.
  for(int k=1;k<=ALPHA;k++){
    const auto&pm=layers[k-1].m;Table tab;tab.build(pm);const auto&cm=layers[k].m;cur.assign(cm.size(),{});
    std::atomic<size_t> next{0};std::array<Accum,THREADS> acc;
    auto worker=[&](int id){
      constexpr size_t BLOCK=2048;
      for(;;){size_t lo=next.fetch_add(BLOCK);if(lo>=cm.size())break;size_t hi=std::min(cm.size(),lo+BLOCK);
        for(size_t z=lo;z<hi;z++){
          Mask s=cm[z],q=s,images[8];for(int op=0;op<8;op++)images[op]=transform(s,op);U256 total{};
          while(q){Mask b=q&-q;q-=b;int v=__builtin_ctzll(b);Mask child=images[0]^mapbit[0][v];for(int op=1;op<8;op++)child=std::min(child,images[op]^mapbit[op][v]);add(total,prev[tab.get(child)]);}
          cur[z]=divsmall(total,k+__builtin_popcountll(neighbourhood(s)&~s));
          if(k>=2){
            const int cuts=n-1;const int orb=layers[k].orb[z];
            int64_t sl[MAXCUT]{}, sh[MAXCUT][MAXCUT]{};
            for(int op=0;op<8;op++){
              int64_t sum[MAXCUT]{}, diagonal[MAXCUT][MAXCUT]{};Mask qq=images[op];
              while(qq){Mask b=qq&-qq;qq-=b;int x=__builtin_ctzll(b)%n;int sv[MAXCUT]{};
                for(int j=0;j<cuts;j++){sv[j]=x<=j?-1:1;sum[j]+=sv[j];}
                for(int j=0;j<cuts;j++)for(int l=0;l<cuts;l++)diagonal[j][l]+=int64_t(sv[j])*sv[l];
              }
              for(int j=0;j<cuts;j++){sl[j]+=int64_t(k-1)*sum[j];for(int l=0;l<cuts;l++)sh[j][l]+=sum[j]*sum[l]-diagonal[j][l];}
            }
            // Orbit average: (orb/8) times the eight transformed values.
            // The product is integral even for stabilised representatives.
            bool neg=(k&1);
            addsigned(acc[id].a2p,acc[id].a2m,cur[z],(neg?-1:1)*int64_t(choose(k,2))*orb);
            for(int j=0;j<cuts;j++){
              assert((sl[j]*orb)%8==0);int64_t q=sl[j]*orb/8;addsigned(acc[id].lp[j],acc[id].lm[j],cur[z],neg?-q:q);
              for(int l=0;l<cuts;l++){assert((sh[j][l]*orb)%8==0);q=sh[j][l]*orb/8;addsigned(acc[id].hp[j][l],acc[id].hm[j][l],cur[z],neg?-q:q);}
            }
          }
        }
      }
    };
    th.clear();for(int i=0;i<THREADS;i++)th.emplace_back(worker,i);for(auto&t:th)t.join();
    for(int i=0;i<THREADS;i++){add(a2p,acc[i].a2p);add(a2m,acc[i].a2m);for(int j=0;j<n-1;j++){add(lp[j],acc[i].lp[j]);add(lm[j],acc[i].lm[j]);for(int l=0;l<n-1;l++){add(hp[j][l],acc[i].hp[j][l]);add(hm[j][l],acc[i].hm[j][l]);}}}
    std::cerr<<"finished k="<<k<<" reps="<<cm.size()<<"\n";
    prev.swap(cur);
  }
  assert(cmp(a2p,a2m)>0);U256 A2=sub(a2p,a2m);const int cuts=n-1;U512 denominator=twice(twice(mul256(A2,A2)));bool nonpositive=true;
  std::cout<<"{\"n\":"<<n<<",\"threads\":"<<THREADS<<",\"raw_independent_sets\":"<<weighted<<",\"representatives_by_size\":[";
  for(int k=0;k<=ALPHA;k++){if(k)std::cout<<',';std::cout<<layers[k].m.size();}std::cout<<"],\"A2\":\""<<str(A2)<<"\",\"cut_marginal_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'\"'<<sstr(lp[j],lm[j])<<'\"';}std::cout<<"],\"ordered_pair_cross_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'[';for(int l=0;l<cuts;l++){if(l)std::cout<<',';std::cout<<'\"'<<sstr(hp[j][l],hm[j][l])<<'\"';}std::cout<<']';}
  std::cout<<"],\"covariance_common_denominator\":\""<<str512(denominator)<<"\",\"covariance_numerators\":[";
  for(int j=0;j<cuts;j++){if(j)std::cout<<',';std::cout<<'[';for(int l=0;l<cuts;l++){if(l)std::cout<<',';SM512 left=product({1,A2},signmag(hp[j][l],hm[j][l]));left.mag=twice(left.mag);SM512 right=product(signmag(lp[j],lm[j]),signmag(lp[l],lm[l]));SM512 cov=subtract_signed(left,right);if(cov.sign>0)nonpositive=false;std::cout<<'\"'<<sstr512(cov)<<'\"';}std::cout<<']';}
  std::cout<<"],\"all_covariance_entries_nonpositive\":"<<(nonpositive?"true":"false")<<"}\n";
}
