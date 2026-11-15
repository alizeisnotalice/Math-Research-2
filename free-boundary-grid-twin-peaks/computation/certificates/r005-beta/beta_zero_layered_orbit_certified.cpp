// Exact beta=0 count on free n-by-n grids, 2<=n<=7.
//
// Unlike beta_zero_orbits.cpp, this program retains only two cardinality
// layers of independent subsets I of U=V\(closed(v) union closed(x)).  The
// stored C, X, E values are <= 49!, hence use packed U224 (84 bytes for all
// three values per state).  Temporary recurrence sums and signed totals use
// U256.  An index-only splitmix64 table maps a predecessor mask to an index:
// keys remain in the predecessor mask array, so the hash adds only 4 bytes
// per slot rather than duplicating 64-bit keys.
//
// The n=7 job is deliberately one process with 12 workers.  It first makes a
// row-mask transfer preflight and rejects an estimated live allocation above
// 8 GiB.  Run `--preflight 7` for the exact layer counts without GMP work.
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using Mask=uint64_t;
using u64=uint64_t;
using u32=uint32_t;
constexpr int THREADS=12;
constexpr Mask EMPTY_MASK=~Mask(0);
constexpr u32 EMPTY_INDEX=std::numeric_limits<u32>::max();

struct U256 { std::array<u64,4> w{}; };
static bool zero(const U256&a){return !(a.w[0]|a.w[1]|a.w[2]|a.w[3]);}
static int cmp(const U256&a,const U256&b){for(int i=3;i>=0;--i)if(a.w[i]!=b.w[i])return a.w[i]<b.w[i]?-1:1;return 0;}
static void add(U256&a,const U256&b){
  u64 carry=0;
  for(int i=0;i<4;i++){u64 z=a.w[i]+b.w[i];bool c1=z<a.w[i];u64 q=z+carry;bool c2=q<z;a.w[i]=q;carry=c1||c2;}
  if(carry)throw std::overflow_error("U256 addition overflow");
}
static U256 sub(U256 a,const U256&b){
  u64 borrow=0;
  for(int i=0;i<4;i++){u64 bi=b.w[i]+borrow;bool c1=bi<b.w[i];u64 old=a.w[i];a.w[i]-=bi;borrow=c1||(old<bi);}
  if(borrow)throw std::overflow_error("negative U256 subtraction"); return a;
}
static U256 mulsmall(const U256&a,u64 q){
  U256 r;__uint128_t carry=0;
  for(int i=0;i<4;i++){__uint128_t z=(__uint128_t)a.w[i]*q+carry;r.w[i]=(u64)z;carry=z>>64;}
  if(carry)throw std::overflow_error("U256 multiply overflow");return r;
}
static U256 divsmall(const U256&a,u64 q){
  U256 r;__uint128_t rem=0;
  for(int i=3;i>=0;i--){__uint128_t z=(rem<<64)|a.w[i];r.w[i]=(u64)(z/q);rem=z%q;}
  if(rem)throw std::runtime_error("nonintegral recurrence division");return r;
}
static std::string dec(U256 a){
  if(zero(a))return "0"; std::string s;
  while(!zero(a)){__uint128_t rem=0;for(int i=3;i>=0;i--){__uint128_t z=(rem<<64)|a.w[i];a.w[i]=(u64)(z/10);rem=z%10;}s.push_back(char('0'+int(rem)));}
  std::reverse(s.begin(),s.end());return s;
}
static U256 factorial(int n){U256 a;a.w[0]=1;for(int i=2;i<=n;i++)a=mulsmall(a,i);return a;}

// Structure-of-arrays avoids the alignment padding a 224-bit struct would
// otherwise incur: 3*64+32 bits exactly per stored value.
struct Layer {
  std::vector<Mask> m;
  std::vector<u64> c0,c1,c2,x0,x1,x2,e0,e1,e2;
  std::vector<u32> c3,x3,e3;
  void resize(size_t n){
    if(n>std::numeric_limits<u32>::max())throw std::runtime_error("too many layer states");
    m.resize(n);c0.resize(n);c1.resize(n);c2.resize(n);c3.resize(n);
    x0.resize(n);x1.resize(n);x2.resize(n);x3.resize(n);
    e0.resize(n);e1.resize(n);e2.resize(n);e3.resize(n);
  }
  size_t size()const{return m.size();}
  U256 getc(size_t i)const{return {{c0[i],c1[i],c2[i],c3[i]}};}
  U256 getx(size_t i)const{return {{x0[i],x1[i],x2[i],x3[i]}};}
  U256 gete(size_t i)const{return {{e0[i],e1[i],e2[i],e3[i]}};}
  static void put(U256 z,std::vector<u64>&a,std::vector<u64>&b,std::vector<u64>&c,std::vector<u32>&d,size_t i){
    if(z.w[3]>>32)throw std::overflow_error("stored value exceeds U224");
    a[i]=z.w[0];b[i]=z.w[1];c[i]=z.w[2];d[i]=(u32)z.w[3];
  }
  void putc(size_t i,const U256&z){put(z,c0,c1,c2,c3,i);}
  void putx(size_t i,const U256&z){put(z,x0,x1,x2,x3,i);}
  void pute(size_t i,const U256&z){put(z,e0,e1,e2,e3,i);}
};

static u64 splitmix64(u64 x){
  x+=0x9e3779b97f4a7c15ULL;
  x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
  x=(x^(x>>27))*0x94d049bb133111ebULL;
  return x^(x>>31);
}

// Slots store only an index.  The actual key comparison uses keys[index].
struct IndexHash {
  std::unique_ptr<std::atomic<u32>[]> slot; const std::vector<Mask>* keys=nullptr; size_t cap=0,mask=0;
  void build(const std::vector<Mask>&a){
    keys=&a;cap=1;while(cap<a.size()*2)cap<<=1; // <=1/2 load, robust for the long n=7 run
    mask=cap-1;slot=std::make_unique<std::atomic<u32>[]>(cap);
    std::atomic<size_t> next{0};std::vector<std::thread> th;
    for(int t=0;t<THREADS;t++)th.emplace_back([&]{for(;;){size_t lo=next.fetch_add(1<<16);if(lo>=cap)break;size_t hi=std::min(cap,lo+(size_t(1)<<16));for(size_t j=lo;j<hi;j++)slot[j].store(0,std::memory_order_relaxed);}});
    for(auto&t:th)t.join();th.clear();next.store(0);
    for(int t=0;t<THREADS;t++)th.emplace_back([&]{for(;;){size_t lo=next.fetch_add(4096);if(lo>=a.size())break;size_t hi=std::min(a.size(),lo+(size_t)4096);for(u32 i=lo;i<hi;i++){size_t h=splitmix64(a[i])&mask;for(;;){u32 empty=0;if(slot[h].compare_exchange_weak(empty,i+1,std::memory_order_relaxed))break;h=(h+1)&mask;}}}});
    for(auto&t:th)t.join();
  }
  u32 get(Mask x)const{size_t h=splitmix64(x)&mask;for(;;){u32 p=slot[h].load(std::memory_order_relaxed);if(!p)throw std::runtime_error("predecessor missing from previous layer");u32 i=p-1;if((*keys)[i]==x)return i;h=(h+1)&mask;}}
};

struct Grid {
  int n,N;std::vector<Mask> adj,closed;
  explicit Grid(int side):n(side),N(side*side),adj(N),closed(N){
    for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=Mask(1)<<u;
    for(int v=0;v<N;v++)closed[v]=adj[v]|(Mask(1)<<v);
  }
  int mapd4(int v,int op)const{int x=v%n,y=v/n;if(op>=4)x=n-1-x;for(int z=0;z<op%4;z++){int q=n-1-y;y=x;x=q;}return y*n+x;}
};
static int pop(Mask x){return __builtin_popcountll(x);} static int bit(Mask x){return __builtin_ctzll(x);}

struct Generator {
  const Grid&g; Mask U; std::vector<std::vector<Mask>> opts; std::vector<std::vector<std::vector<u64>>> memo;
  struct Task{Mask r0,r1;int used,need;u64 count,offset;};
  Generator(const Grid&G,Mask u):g(G),U(u),opts(g.n),memo(g.n+1){
    for(int y=0;y<g.n;y++){Mask allow=0;for(int x=0;x<g.n;x++){int v=y*g.n+x;if(U&(Mask(1)<<v))allow|=Mask(1)<<x;}
      for(Mask r=0;r<(Mask(1)<<g.n);r++)if(!(r&~allow)&&!(r&(r<<1)))opts[y].push_back(r);}
    for(int y=0;y<=g.n;y++)memo[y].assign(Mask(1)<<g.n,std::vector<u64>(g.N+1,~u64(0)));
  }
  u64 ways(int y,Mask last,int need){
    if(need<0)return 0;if(y==g.n)return need==0;
    u64& z=memo[y][last][need];if(z!=~u64(0))return z;z=0;
    for(Mask r:opts[y])if(!(r&last)){int q=pop(r);if(q<=need)z+=ways(y+1,r,need-q);}return z;
  }
  std::vector<u64> counts(){std::vector<u64>a(g.N+1);for(int k=0;k<=g.N;k++)a[k]=ways(0,0,k);while(a.size()>1&&!a.back())a.pop_back();return a;}
  void emit(int y,Mask last,int need,Mask chosen,Mask* out,u64& at){
    if(y==g.n){if(!need)out[at++]=chosen;return;}
    for(Mask r:opts[y])if(!(r&last)){int q=pop(r);if(q<=need&&ways(y+1,r,need-q))emit(y+1,r,need-q,chosen|(r<<(y*g.n)),out,at);}
  }
  void make_layer(int k,Layer&out){
    std::vector<Task> tasks;
    // Split on first two rows: enough work packets for 12 workers, and it
    // gives exact offsets so no duplicate per-thread mask buffers are needed.
    for(Mask r0:opts[0])for(Mask r1:opts[1])if(!(r0&r1)){int q=pop(r0)+pop(r1);if(q<=k){u64 c=ways(2,r1,k-q);if(c)tasks.push_back({r0,r1,q,k-q,c,0});}}
    u64 total=0;for(auto&t:tasks){t.offset=total;total+=t.count;}out.resize(total);
    std::atomic<size_t> next{0};std::vector<std::thread> th;
    for(int z=0;z<THREADS;z++)th.emplace_back([&]{for(;;){size_t i=next.fetch_add(1);if(i>=tasks.size())break;const auto&t=tasks[i];u64 at=t.offset;emit(2,t.r1,t.need,t.r0|(t.r1<<g.n),out.m.data(),at);if(at!=t.offset+t.count)throw std::runtime_error("layer enumeration count mismatch");}});
    for(auto&t:th)t.join();
  }
};

struct SignedAcc{U256 plus,minus;void addterm(const U256&z,bool positive){if(positive)add(plus,z);else add(minus,z);}U256 value()const{return sub(plus,minus);}};
struct EdgeAnswer{U256 e2,vhigh;std::vector<u64> layers;};

static EdgeAnswer edge_answer(const Grid&g,int v,int x,bool preflight_only){
  Mask T=g.closed[v]|g.closed[x], full=(Mask(1)<<g.N)-1,U=full^T;
  Generator gen(g,U);auto counts=gen.counts();
  if(preflight_only)return {{},{},counts};
  // Exact live-byte upper estimate: two packed layers (92 bytes/state), a
  // half-loaded index hash (4 bytes/slot), plus 512 MiB allocator/thread headroom.
  size_t worst=0;for(size_t k=1;k<counts.size();k++){size_t cap=1;while(cap<2*counts[k-1])cap<<=1;size_t bytes=92*(counts[k-1]+counts[k])+4*cap+(size_t(512)<<20);worst=std::max(worst,bytes);} 
  constexpr size_t LIMIT=size_t(8)<<30;
  if(worst>LIMIT)throw std::runtime_error("preflight exceeds 8 GiB live allocation estimate");
  std::cerr<<"edge ["<<v<<","<<x<<"] states=";for(size_t k=0;k<counts.size();k++)if(counts[k])std::cerr<<k<<":"<<counts[k]<<" ";std::cerr<<" estimated_peak_GiB="<<std::fixed<<std::setprecision(3)<<(double(worst)/(1ULL<<30))<<"\n";
  U256 fact=factorial(g.N);Layer prev;SignedAcc ae,av;
  for(size_t k=0;k<counts.size();k++){
    Layer cur;gen.make_layer((int)k,cur);if(cur.size()!=counts[k])throw std::runtime_error("preflight/layer size disagreement");
    IndexHash hash;if(k)hash.build(prev.m);
    struct Local{U256 ep,em,vp,vm;};std::array<Local,THREADS> local{};std::atomic<size_t> next{0};std::vector<std::thread> th;
    for(int tid=0;tid<THREADS;tid++)th.emplace_back([&,tid]{
      constexpr size_t BLOCK=4096;
      for(;;){size_t lo=next.fetch_add(BLOCK);if(lo>=cur.size())break;size_t hi=std::min(cur.size(),lo+BLOCK);
        for(size_t z=lo;z<hi;z++){Mask I=cur.m[z],q=I,cl=0;U256 sc{},sx{},se{};
          while(q){Mask b=q&-q;q-=b;int i=bit(b);cl|=g.closed[i];u32 p=hash.get(I^b);add(sc,prev.getc(p));add(sx,prev.getx(p));add(se,prev.gete(p));}
          U256 C;if(!I)C=fact;else C=divsmall(sc,pop(cl));
          int a=pop(cl|T);U256 X=C;add(X,sx);X=divsmall(X,a-1);U256 E=X;add(E,se);E=divsmall(E,a);
          cur.putc(z,C);cur.putx(z,X);cur.pute(z,E);
          if(I){U256 VX=divsmall(X,g.N),qE=mulsmall(E,k),qX=mulsmall(VX,k);if(k&1){add(local[tid].ep,qE);add(local[tid].vp,qX);}else{add(local[tid].em,qE);add(local[tid].vm,qX);}}
        }
      }
    });
    for(auto&t:th)t.join();
    for(auto&z:local){add(ae.plus,z.ep);add(ae.minus,z.em);add(av.plus,z.vp);add(av.minus,z.vm);} prev=std::move(cur);
  }
  U256 e2=ae.value(), vh=av.value();if(cmp(e2,vh)<0)throw std::runtime_error("negative beta0 edge count");return {e2,vh,counts};
}

int main(int argc,char**argv){
 try{
  bool pre=false;int n=0,chosen_orbit=-1;
  if(argc==3&&std::string(argv[1])=="--preflight"){pre=true;n=std::stoi(argv[2]);}
  else if(argc==2)n=std::stoi(argv[1]);
  else if(argc==4&&std::string(argv[2])=="--orbit"){n=std::stoi(argv[1]);chosen_orbit=std::stoi(argv[3]);}
  else return 2;
  if(n<2||n>7)return 2;Grid g(n);std::map<std::pair<int,int>,std::vector<std::pair<int,int>>> bucket;
  for(int v=0;v<g.N;v++)for(Mask q=g.adj[v];q;q&=q-1){int x=bit(q);if(v>=x)continue;std::pair<int,int> key={g.N,g.N};for(int op=0;op<8;op++){int a=g.mapd4(v,op),b=g.mapd4(x,op);if(a>b)std::swap(a,b);key=std::min(key,{a,b});}bucket[key].push_back({v,x});}
  size_t undirected=0;for(const auto&[_,o]:bucket)undirected+=o.size();if(undirected!=size_t(2*n*(n-1)))throw std::runtime_error("D4 undirected orbit coverage mismatch");
  if(chosen_orbit>=(int)bucket.size())return 2;
  U256 total{};size_t selected_undirected=0,orbit_number=0;std::cout<<"{\"n\":"<<n<<",\"edge_orbits\":[";bool first=true;
  for(const auto&[key,o]:bucket){if(chosen_orbit>=0&&int(orbit_number++)!=chosen_orbit)continue;auto ans=edge_answer(g,o.front().first,o.front().second,pre);if(!first)std::cout<<',';first=false;std::cout<<"{\"representative\":["<<o.front().first<<','<<o.front().second<<"],\"multiplicity\":"<<o.size()<<",\"layers\":[";for(size_t i=0;i<ans.layers.size();i++){if(i)std::cout<<',';std::cout<<ans.layers[i];}std::cout<<']';selected_undirected+=o.size();if(!pre){U256 low=sub(ans.e2,ans.vhigh);add(total,mulsmall(low,2*o.size()));std::cout<<",\"E_K2\":\""<<dec(ans.e2)<<"\",\"v_high\":\""<<dec(ans.vhigh)<<"\",\"beta0_one_direction\":\""<<dec(low)<<'"';}std::cout<<'}';}
  std::cout<<"],\"undirected_edge_coverage\":"<<selected_undirected<<",\"directed_edge_coverage\":"<<2*selected_undirected<<",\"global_directed_edge_coverage\":"<<2*undirected;if(chosen_orbit>=0)std::cout<<",\"selected_orbit\":"<<chosen_orbit;if(!pre)std::cout<<",\"beta0_total\":\""<<dec(total)<<'"';std::cout<<"}\n";return 0;
 }catch(const std::exception&e){std::cerr<<"beta_zero_layered: "<<e.what()<<"\n";return 1;}
}
