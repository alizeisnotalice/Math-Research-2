// Exact n=7 peak-set inclusion-exclusion, reduced by the D4 action.
// Every stored mask is the least member of its D4 orbit.  C(I)=49! P(I all
// peaks) is computed cardinality layer by layer from C(I\{v}).  The program
// uses one multithreaded process; its hard memory budget is well below 10 GiB.
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
constexpr int n=7, N=49, ALPHA=25, THREADS=12;
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

static std::array<Mask,N> closedv;
static std::array<std::array<int,N>,N> dist;
static std::array<Mask,N> mapbit[8];
static std::vector<Mask> rows;
static uint64_t raw_by_k[ALPHA+1]={
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
static std::array<Layer,ALPHA+1> layers;
static std::atomic<int> next_task{0};
static std::vector<std::pair<Mask,Mask>> starts;
static uint64_t hash64(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;return x;}

static void enum_tail(int row, Mask previous, Mask set, std::array<std::vector<Mask>,ALPHA+1>&localm,
                      std::array<std::vector<uint8_t>,ALPHA+1>&localo){
  if(row==n){uint8_t o;if(canon(set,&o)==set){int k=__builtin_popcountll(set);localm[k].push_back(set);localo[k].push_back(o);}return;}
  for(Mask r:rows)if(!(r&previous))enum_tail(row+1,r,set|(r<<(row*n)),localm,localo);
}
static void enumerate_worker(){
  std::array<std::vector<Mask>,ALPHA+1> lm;
  std::array<std::vector<uint8_t>,ALPHA+1> lo;
  for(;;){int t=next_task.fetch_add(1);if(t>=(int)starts.size())break;auto [a,b]=starts[t];enum_tail(2,b,a|(b<<n),lm,lo);}
  static std::mutex mu;
  std::lock_guard<std::mutex> g(mu);
  for(int k=0;k<=ALPHA;k++){
    layers[k].m.insert(layers[k].m.end(),lm[k].begin(),lm[k].end());
    layers[k].orb.insert(layers[k].orb.end(),lo[k].begin(),lo[k].end());
  }
}

// The immutable previous-layer masks carry the hash keys.  Slots store index+1,
// so a 32-bit atomic slot is enough and construction can be parallelized.
struct Table{
  std::vector<std::atomic<uint32_t>> slot; const std::vector<Mask>* keys=nullptr; uint64_t mask=0;
  double clear_seconds=0,insert_seconds=0;
  void build(const std::vector<Mask>&a){
    keys=&a;uint64_t cap=1;while(cap<a.size()*2)cap<<=1;slot=std::vector<std::atomic<uint32_t>>(cap);mask=cap-1;
    auto parallel=[&](size_t total,auto body){std::atomic<size_t> next{0};std::vector<std::thread> ts;for(int q=0;q<THREADS;q++)ts.emplace_back([&]{for(;;){size_t lo=next.fetch_add(1<<16);if(lo>=total)break;body(lo,std::min(total,lo+(1<<16)));}});for(auto&t:ts)t.join();};
    auto start=std::chrono::steady_clock::now();
    parallel(cap,[&](size_t lo,size_t hi){for(size_t i=lo;i<hi;i++)slot[i].store(0,std::memory_order_relaxed);});
    clear_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();start=std::chrono::steady_clock::now();
    parallel(a.size(),[&](size_t lo,size_t hi){for(size_t i=lo;i<hi;i++){uint64_t h=hash64(a[i])&mask;uint32_t wanted=uint32_t(i+1);for(;;){uint32_t old=slot[h].load(std::memory_order_relaxed);if(!old){if(slot[h].compare_exchange_weak(old,wanted,std::memory_order_relaxed,std::memory_order_relaxed))break;}else {assert((*keys)[old-1]!=a[i]);h=(h+1)&mask;}}}});
    insert_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  }
  uint32_t get(Mask x)const{uint64_t h=hash64(x)&mask;for(;;){uint32_t p=slot[h].load(std::memory_order_relaxed);assert(p);if((*keys)[p-1]==x)return p-1;h=(h+1)&mask;}}
};

static Mask neighbourhood(Mask s){Mask z=0;while(s){Mask b=s&-s;s-=b;z|=closedv[__builtin_ctzll(b)];}return z;}
static uint64_t choose(int a,int b){if(b<0||b>a)return 0;b=std::min(b,a-b);uint64_t x=1;for(int i=1;i<=b;i++)x=x*(a-b+i)/i;return x;}
struct Accum{std::array<U256,ALPHA+1> t{};std::array<U256,2*N-1> d{};U256 sq{};};

int main(){
  for(int y=0;y<n;y++)for(int x=0;x<n;x++){
    int v=y*n+x;closedv[v]=Mask(1)<<v;
    for(auto [X,Y]:{std::pair{x-1,y},{x+1,y},{x,y-1},{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)closedv[v]|=Mask(1)<<(Y*n+X);
    for(int w=0;w<N;w++)dist[v][w]=std::abs(x-w%n)+std::abs(y-w/n);
    for(int op=0;op<8;op++){int X=x,Y=y;if(op>=4)X=n-1-X;for(int q=0;q<op%4;q++){int T=n-1-Y;Y=X;X=T;}mapbit[op][v]=Mask(1)<<(Y*n+X);}
  }
  for(Mask r=0;r<(1<<n);r++)if(!(r&(r<<1)))rows.push_back(r);
  for(Mask a:rows)for(Mask b:rows)if(!(a&b))starts.push_back({a,b});
  std::cerr<<"enumerating "<<starts.size()<<" two-row tasks with "<<THREADS<<" threads\n";
  std::vector<std::thread> th;for(int i=0;i<THREADS;i++)th.emplace_back(enumerate_worker);for(auto&t:th)t.join();
  uint64_t weighted=0;
  for(int k=0;k<=ALPHA;k++){uint64_t got=0;for(auto o:layers[k].orb)got+=o;weighted+=got;assert(got==raw_by_k[k]);std::cerr<<"k="<<k<<" reps="<<layers[k].m.size()<<" weighted="<<got<<"\n";}
  assert(weighted==1280128950ULL);
  std::vector<U256> prev(1),cur;
  U256 fact;fact.w[0]=1;for(int i=2;i<=N;i++)fact=mulsmall(fact,i);prev[0]=fact;
  std::array<Mask,4> probes{0,Mask(1),Mask(1)|(Mask(1)<<2),0};
  for(int y=0;y<n;y+=2)for(int x=0;x<n;x+=2)probes[3]|=Mask(1)<<(y*n+x);
  std::array<U256,4> probeC{};probeC[0]=fact;
  std::array<U256,ALPHA+1> TkP{},TkM{};std::array<U256,2*N-1> Dp{},Dm{};U256 Sp{},Sm{};
  double hash_clear_total=0,hash_insert_total=0,recur_total=0;
  // k=0 contributes neither peaks nor pairs.  All other layers are independent.
  for(int k=1;k<=ALPHA;k++){
    const auto&pm=layers[k-1].m;Table tab;tab.build(pm);hash_clear_total+=tab.clear_seconds;hash_insert_total+=tab.insert_seconds;const auto&cm=layers[k].m;cur.assign(cm.size(),{});
    auto recur_start=std::chrono::steady_clock::now();
    std::atomic<size_t> next{0};std::array<Accum,THREADS> acc;
    auto worker=[&](int id){
      constexpr size_t BLOCK=2048;
      for(;;){size_t lo=next.fetch_add(BLOCK);if(lo>=cm.size())break;size_t hi=std::min(cm.size(),lo+BLOCK);
        for(size_t z=lo;z<hi;z++){
          Mask s=cm[z],q=s,images[8];for(int op=0;op<8;op++)images[op]=transform(s,op);U256 total{};
          while(q){Mask b=q&-q;q-=b;int v=__builtin_ctzll(b);Mask child=images[0]^mapbit[0][v];for(int op=1;op<8;op++)child=std::min(child,images[op]^mapbit[op][v]);add(total,prev[tab.get(child)]);}
          cur[z]=divsmall(total,k+__builtin_popcountll(neighbourhood(s)&~s));
          U256 w=mulsmall(cur[z],layers[k].orb[z]); add(acc[id].t[k],w);
          if(k>=2){int vs[ALPHA],nv=0;Mask qq=s;while(qq){Mask b=qq&-qq;qq-=b;vs[nv++]=__builtin_ctzll(b);}uint64_t ss=0;
            for(int i=0;i<nv;i++)for(int j=i+1;j<nv;j++){int dd=dist[vs[i]][vs[j]];add(acc[id].d[dd],w);ss+=uint64_t((vs[i]%n-vs[j]%n)*(vs[i]%n-vs[j]%n)+(vs[i]/n-vs[j]/n)*(vs[i]/n-vs[j]/n));}
            add(acc[id].sq,mulsmall(w,ss));}
        }
      }
    };
    th.clear();for(int i=0;i<THREADS;i++)th.emplace_back(worker,i);for(auto&t:th)t.join();recur_total+=std::chrono::duration<double>(std::chrono::steady_clock::now()-recur_start).count();
    for(int i=0;i<THREADS;i++){add(TkP[k],acc[i].t[k]);for(int d=0;d<2*N-1;d++){if(k&1)add(Dm[d],acc[i].d[d]);else add(Dp[d],acc[i].d[d]);}if(k&1)add(Sm,acc[i].sq);else add(Sp,acc[i].sq);}
    for(int p=1;p<4;p++)if(__builtin_popcountll(probes[p])==k){Table lookup;lookup.build(cm);probeC[p]=cur[lookup.get(canon(probes[p]))];}
    // Histogram signs are applied after all layer totals are known.
    std::cerr<<"finished k="<<k<<" reps="<<cm.size()<<" hash_clear_s="<<tab.clear_seconds<<" hash_insert_s="<<tab.insert_seconds<<" recur_s="<<recur_total<<"\n";
    prev.swap(cur);
  }
  std::array<U256,ALPHA+1> hp{},hm{};
  for(int r=1;r<=ALPHA;r++)for(int k=r;k<=ALPHA;k++){U256 q=mulsmall(TkP[k],choose(k,r));if((k-r)&1)add(hm[r],q);else add(hp[r],q);}
  std::cerr<<"stage_totals hash_clear_s="<<hash_clear_total<<" hash_insert_s="<<hash_insert_total<<" recurrence_s="<<recur_total<<"\n";
  std::cout<<"{\"n\":7,\"raw_independent_sets\":1280128950,\"representatives_by_size\":[";
  for(int k=0;k<=ALPHA;k++){if(k)std::cout<<',';std::cout<<layers[k].m.size();}std::cout<<"],\"peak_histogram\":{";
  for(int k=1;k<=ALPHA;k++){if(k>1)std::cout<<',';std::cout<<'\"'<<k<<"\":\""<<sstr(hp[k],hm[k])<<"\"";}std::cout<<"},\"A1\":\""<<sstr(hp[1],hm[1])<<"\",\"A2\":\""<<sstr(hp[2],hm[2])<<"\",\"distance_counts\":{";
  bool first=true;for(int d=0;d<=12;d++)if(cmp(Dp[d],Dm[d])){if(!first)std::cout<<',';first=false;std::cout<<'\"'<<d<<"\":\""<<sstr(Dp[d],Dm[d])<<"\"";}std::cout<<"},\"squared_euclidean_total\":\""<<sstr(Sp,Sm)<<"\",\"selected_joint_peak_counts\":{\"empty\":\""<<str(probeC[0])<<"\",\"corner\":\""<<str(probeC[1])<<"\",\"two_same_row\":\""<<str(probeC[2])<<"\",\"even_even_4x4\":\""<<str(probeC[3])<<"\"}}\n";
}
