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
constexpr int MAXN=49, MAXALPHA=25, THREADS=4;
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

static std::array<Mask,MAXN> closedv;
static std::array<std::array<int,MAXN>,MAXN> dist;
static std::array<Mask,MAXN> mapbit[8];
static std::array<uint8_t,MAXN> vclass, classsize;
static std::array<Mask,MAXN> classmask;
static int classes=0;
constexpr int MAXSEEDS=8;
static std::vector<Mask> seeds, bans;
static std::vector<std::string> seed_names;
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
struct Accum{
  std::array<U256,MAXALPHA+1> t{};
  std::array<U256,2*MAXN-1> d{};
  std::array<U256,MAXN> rp{},rm{};
  U256 sq{},fp{},fm{},gp{},gm{};
  std::array<U256,MAXSEEDS> seedp{},seedm{};
};

int main(int argc,char**argv){
  if(argc!=2)return 2; n=std::atoi(argv[1]); if(n<3||n>7)return 2;
  N=n*n; ALPHA=(N+1)/2;
  for(int y=0;y<n;y++)for(int x=0;x<n;x++){
    int v=y*n+x;closedv[v]=Mask(1)<<v;
    for(auto [X,Y]:{std::pair{x-1,y},{x+1,y},{x,y-1},{x,y+1}})if(0<=X&&X<n&&0<=Y&&Y<n)closedv[v]|=Mask(1)<<(Y*n+X);
    for(int w=0;w<N;w++)dist[v][w]=std::abs(x-w%n)+std::abs(y-w/n);
    for(int op=0;op<8;op++){int X=x,Y=y;if(op>=4)X=n-1-X;for(int q=0;q<op%4;q++){int T=n-1-Y;Y=X;X=T;}mapbit[op][v]=Mask(1)<<(Y*n+X);}
  }
  vclass.fill(255);
  for(int v=0;v<N;v++)if(vclass[v]==255){
    Mask orbit=0;
    for(int op=0;op<8;op++)orbit|=mapbit[op][v];
    int c=classes++, size=__builtin_popcountll(orbit);
    classmask[c]=orbit; classsize[c]=size;
    for(Mask q=orbit;q;){Mask b=q&-q;q-=b;vclass[__builtin_ctzll(b)]=c;}
  }
  for(int v=0;v<N;v++)assert(vclass[v]<classes);
  auto seed=[&](Mask m,std::string name){seeds.push_back(m);bans.push_back(neighbourhood(m));seed_names.push_back(name);};
  seed(Mask(1),"corner-root");seed(Mask(1)<<(n+1),"inner-root");
  for(int a=1;a<=std::min(3,n-2);a++){
    Mask m=1;for(int y=1;y<=a;y++)for(int x=1;x<=a;x++)m|=Mask(1)<<(y*n+x);
    seed(m,"block"+std::to_string(a)+"-dormant");seed(m|2,"block"+std::to_string(a)+"-bridge");
  }
  assert(seeds.size()<=MAXSEEDS);
  for(Mask r=0;r<(1<<n);r++)if(!(r&(r<<1)))rows.push_back(r);
  for(Mask a:rows)for(Mask b:rows)if(!(a&b))starts.push_back({a,b});
  std::cerr<<"enumerating "<<starts.size()<<" two-row tasks with "<<THREADS<<" threads\n";
  std::vector<std::thread> th;for(int i=0;i<THREADS;i++)th.emplace_back(enumerate_worker);for(auto&t:th)t.join();
  uint64_t weighted=0;
  for(int k=0;k<=ALPHA;k++){uint64_t got=0;for(auto o:layers[k].orb)got+=o;weighted+=got;if(n==7)assert(got==raw_by_k[k]);std::cerr<<"k="<<k<<" reps="<<layers[k].m.size()<<" weighted="<<got<<"\n";}
  if(n==7)assert(weighted==1280128950ULL);
  std::vector<U256> prev(1),cur;
  U256 fact;fact.w[0]=1;for(int i=2;i<=N;i++)fact=mulsmall(fact,i);prev[0]=fact;
  std::array<U256,MAXALPHA+1> TkP{},TkM{};
  std::array<U256,2*MAXN-1> Dp{},Dm{};
  std::array<U256,MAXN> RootP{},RootM{};
  U256 Sp{},Sm{},Fp{},Fm{},Gp{},Gm{};
  std::array<U256,MAXSEEDS> SeedP{},SeedM{};for(size_t j=0;j<seeds.size();j++)SeedP[j]=fact;
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
          U256 w=mulsmall(cur[z],layers[k].orb[z]); add(acc[id].t[k],w);
          // Every distinct orbit member occurs 8/orbit_size times among
          // images[0..7]. Multiply BEFORE division to retain stabilizers.
          for(size_t j=0;j<seeds.size();j++){
            unsigned hits=0;for(int op=0;op<8;op++)hits+=!(images[op]&bans[j]);
            unsigned mult=hits*layers[k].orb[z];assert(mult%8==0);mult/=8;
            if(mult){U256 term=mulsmall(cur[z],mult);if(k&1)add(acc[id].seedm[j],term);else add(acc[id].seedp[j],term);}
          }

          // For one fixed root in a D4 vertex orbit J, precisely
          // orbit(s)*|s cap J|/|J| members of the independent-set orbit
          // contain it.  This is an integer by transitivity of D4 on J.
          for(int c=0;c<classes;c++){
            int hits=__builtin_popcountll(s&classmask[c]);
            if(!hits)continue;
            uint64_t q=uint64_t(layers[k].orb[z])*hits;
            assert(q%classsize[c]==0); q/=classsize[c];
            U256 term=mulsmall(cur[z],q);
            for(Mask u=classmask[c];u;){Mask b=u&-u;u-=b;int v=__builtin_ctzll(b);
              if(k&1)add(acc[id].rp[v],term);else add(acc[id].rm[v],term);
            }
          }
          if(k>=2){int vs[MAXALPHA],nv=0;Mask qq=s;while(qq){Mask b=qq&-qq;qq-=b;vs[nv++]=__builtin_ctzll(b);}uint64_t ss=0;
            for(int i=0;i<nv;i++)for(int j=i+1;j<nv;j++){int dd=dist[vs[i]][vs[j]];add(acc[id].d[dd],w);ss+=uint64_t((vs[i]%n-vs[j]%n)*(vs[i]%n-vs[j]%n)+(vs[i]/n-vs[j]/n)*(vs[i]/n-vs[j]/n));}
            add(acc[id].sq,mulsmall(w,ss));}
          if(k>=2){
            int64_t sx=0,sy=0,norm=0;Mask qq=s;
            while(qq){Mask b=qq&-qq;qq-=b;int v=__builtin_ctzll(b);int64_t tx=2*(v%n)-(n-1),ty=2*(v/n)-(n-1);sx+=tx;sy+=ty;norm+=tx*tx+ty*ty;}
            int64_t f=(k-1)*norm,g=(sx*sx+sy*sy-norm)/2;assert((sx*sx+sy*sy-norm)%2==0);
            bool neg=(k&1);if(neg)add(acc[id].fm,mulsmall(w,f));else add(acc[id].fp,mulsmall(w,f));
            if(g<0){neg=!neg;g=-g;}if(neg)add(acc[id].gm,mulsmall(w,g));else add(acc[id].gp,mulsmall(w,g));
          }
        }
      }
    };
    th.clear();for(int i=0;i<THREADS;i++)th.emplace_back(worker,i);for(auto&t:th)t.join();
    for(int i=0;i<THREADS;i++){
      add(TkP[k],acc[i].t[k]);
      for(size_t j=0;j<seeds.size();j++){add(SeedP[j],acc[i].seedp[j]);add(SeedM[j],acc[i].seedm[j]);}
      for(int v=0;v<N;v++){add(RootP[v],acc[i].rp[v]);add(RootM[v],acc[i].rm[v]);}
      for(int d=0;d<2*N-1;d++){if(k&1)add(Dm[d],acc[i].d[d]);else add(Dp[d],acc[i].d[d]);}
      if(k&1)add(Sm,acc[i].sq);else add(Sp,acc[i].sq);
      add(Fp,acc[i].fp);add(Fm,acc[i].fm);add(Gp,acc[i].gp);add(Gm,acc[i].gm);
    }
    // Histogram signs are applied after all layer totals are known.
    std::cerr<<"finished k="<<k<<" reps="<<cm.size()<<"\n";
    prev.swap(cur);
  }
  std::array<U256,MAXALPHA+1> hp{},hm{};
  for(int r=1;r<=ALPHA;r++)for(int k=r;k<=ALPHA;k++){U256 q=mulsmall(TkP[k],choose(k,r));if((k-r)&1)add(hm[r],q);else add(hp[r],q);}
  U256 rootpos{},rootneg{};
  for(int v=0;v<N;v++){add(rootpos,RootP[v]);add(rootneg,RootM[v]);}
  assert(sstr(rootpos,rootneg)==sstr(hp[1],hm[1]));
  for(int c=0;c<classes;c++){Mask q=classmask[c];int v=__builtin_ctzll(q);for(q&=q-1;q;q&=q-1)assert(sstr(RootP[v],RootM[v])==sstr(RootP[__builtin_ctzll(q)],RootM[__builtin_ctzll(q)]));}
  std::cout<<"{\"fixed_state_completions\":[";
  for(size_t j=0;j<seeds.size();j++){
    assert(cmp(SeedP[j],SeedM[j])>=0);U256 h=sub(SeedP[j],SeedM[j]);
    unsigned sz=__builtin_popcountll(seeds[j]);for(unsigned t=0;t<sz;t++)h=divsmall(h,N-t);
    if(j==0)assert(str(h)==sstr(RootP[0],RootM[0]));
    if(j==1)assert(str(h)==sstr(RootP[n+1],RootM[n+1]));
    if(j)std::cout<<',';
    std::cout<<"{\"name\":\""<<seed_names[j]<<"\",\"S\":"<<seeds[j]<<",\"size\":"<<sz<<",\"ban\":"<<bans[j]<<",\"completion_count\":\""<<str(h)<<"\"}";
  }
  std::cout<<"],";
  std::cout<<"\"n\":"<<n<<",\"threads\":"<<THREADS<<",\"raw_independent_sets\":"<<weighted<<",\"representatives_by_size\":[";
  for(int k=0;k<=ALPHA;k++){if(k)std::cout<<',';std::cout<<layers[k].m.size();}std::cout<<"],\"peak_histogram\":{";
  for(int k=1;k<=ALPHA;k++){if(k>1)std::cout<<',';std::cout<<'\"'<<k<<"\":\""<<sstr(hp[k],hm[k])<<"\"";}
  std::cout<<"},\"A1\":\""<<sstr(hp[1],hm[1])<<"\",\"root_counts\":[";
  for(int v=0;v<N;v++){if(v)std::cout<<',';std::cout<<'\"'<<sstr(RootP[v],RootM[v])<<"\"";}
  std::cout<<"],\"A2\":\""<<sstr(hp[2],hm[2])<<"\",\"distance_counts\":{";
  bool first=true;for(int d=0;d<=2*(n-1);d++)if(cmp(Dp[d],Dm[d])){if(!first)std::cout<<',';first=false;std::cout<<'\"'<<d<<"\":\""<<sstr(Dp[d],Dm[d])<<"\"";}std::cout<<"},\"squared_euclidean_total\":\""<<sstr(Sp,Sm)<<"\",\"firstsum\":\""<<sstr(Fp,Fm)<<"\",\"secondsum\":\""<<sstr(Gp,Gm)<<"\"}\n";
}
