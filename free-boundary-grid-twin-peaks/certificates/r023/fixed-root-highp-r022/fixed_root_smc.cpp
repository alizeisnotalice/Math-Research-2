// Fixed-root, ordered-prefix resample-move SMC for ONE peak only.
// Z estimates h_r: the number of labelled one-peak orderings whose first
// (highest-labelled peak) vertex is root r.  It is an integer-count normalizer,
// not h_r/N! and not a probability.  This r021 copy has no fibre observable.
#include <array>
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>
using U=uint64_t; using H=__uint128_t;
struct W { std::array<unsigned char,121> a{}; };
struct R { double z=1,miness=0,acc=0,props=0; int dead=0; };
struct G {
  int n,N; std::array<H,121> adj{};
  explicit G(int q):n(q),N(q*q) { for(int y=0;y<n;y++) for(int x=0;x<n;x++) {
    int k=y*n+x; for(auto [dx,dy]:{std::pair{1,0},std::pair{-1,0},std::pair{0,1},std::pair{0,-1}}) {
      int X=x+dx,Y=y+dy; if(X>=0&&X<n&&Y>=0&&Y<n) adj[k]|=H(1)<<(Y*n+X);
    }
  }};
};
static int uni(std::mt19937_64& g,int m) { assert(m>0); U cut=(U(0)-U(m))%U(m); for(;;) { U z=g(); if(z>=cut) return int(z%U(m)); } }
static bool valid(const W&w,int r,const G&g) {
  H seen=0; int births=0;
  for(int t=0;t<r;t++) { int v=w.a[t]; H b=H(1)<<v; if(seen&b) return false; if(!(g.adj[v]&seen)) { births++; if(t!=0) return false; } seen|=b; }
  return births==1;
}
static int allowed(const W&w,int r,const G&g,unsigned char*out) {
  H seen=0,bd=0; for(int t=0;t<r;t++) { int v=w.a[t]; seen|=H(1)<<v; bd|=g.adj[v]; }
  H all=(H(1)<<g.N)-1,m=(r==0?all:(bd&~seen)); int k=0;
  for(int v=0;v<g.N;v++) if(m&(H(1)<<v)) out[k++]=v;
  return k;
}
static int one_peak(const W&w,const G&g) {
  std::array<int,121> rank{}; for(int t=0;t<g.N;t++) rank[w.a[t]]=t; int total=0;
  for(int v=0;v<g.N;v++) { bool peak=true; for(int z=0;z<g.N;z++) if((g.adj[v]&(H(1)<<z))&&rank[z]<rank[v]) peak=false; if(peak) total++; }
  return total;
}
// The unrestricted proposal is symmetric on ordered distinct prefixes.  We
// retain it, and reject every endpoint whose first entry is not fixed_root;
// hence the endpoint-rejected kernel is symmetric on the fixed-root support.
static W mh(W w,int r,int fixed_root,int L,const G&g,std::mt19937_64&rg,R&z) {
  for(int q=0;q<L;q++) { z.props++; int i=uni(rg,r), v=uni(rg,g.N), j=-1;
    for(int t=0;t<r;t++) if(w.a[t]==v) { j=t; break; }
    W y=w; if(j>=0) std::swap(y.a[i],y.a[j]); else y.a[i]=v;
    if(y.a[0]==fixed_root && valid(y,r,g)) { w=y; z.acc++; }
  }
  return w;
}
static R run(const G&g,int fixed_root,int P,int L,U seed) {
  assert(0<=fixed_root&&fixed_root<g.N); std::mt19937_64 rg(seed); std::vector<W>x(P);
  for(auto&w:x) w.a[0]=unsigned(fixed_root); R o; unsigned char al[121];
  for(int r=2;r<=g.N;r++) {
    std::vector<int>a(P); double sw=0,s2=0;
    for(int i=0;i<P;i++) { assert(x[i].a[0]==fixed_root && valid(x[i],r-1,g)); a[i]=allowed(x[i],r-1,g,al); sw+=a[i]; s2+=double(a[i])*a[i]; if(!a[i])o.dead++; }
    o.miness=(r==2?sw*sw/s2:std::min(o.miness,sw*sw/s2)); double mz=sw/P; o.z*=mz;
    if(!sw) { o.z=0; return o; }
    if(r==g.N) { for(int i=0;i<P;i++) if(a[i]) { int k=allowed(x[i],r-1,g,al); assert(k==a[i]); W y=x[i]; y.a[r-1]=al[uni(rg,k)]; assert(y.a[0]==fixed_root && valid(y,r,g) && one_peak(y,g)==1); } return o; }
    std::vector<U> c(P); U total=0; for(int i=0;i<P;i++) { assert(UINT64_MAX-total>=U(a[i])); total+=U(a[i]); c[i]=total; }
    std::vector<W> nx; nx.reserve(P);
    for(int q=0;q<P;q++) { U u=U(uni(rg,int(total))); int j=int(std::upper_bound(c.begin(),c.end(),u)-c.begin()); if(j>=P)j=P-1; int k=allowed(x[j],r-1,g,al); assert(k==a[j]); W y=x[j]; y.a[r-1]=al[uni(rg,k)]; y=mh(y,r,fixed_root,L,g,rg,o); assert(y.a[0]==fixed_root && valid(y,r,g)); nx.push_back(y); }
    x.swap(nx);
  }
  return o;
}
// Independent exact subset recurrence: number of growth orders with only root
// r allowed to be born.  It is used only for n<=4, where all values fit U.
static U exact_root(const G&g,int root) {
  assert(g.N<=16); size_t lim=size_t(1)<<g.N; std::vector<U>dp(lim); dp[0]=1; H full=(H(1)<<g.N)-1;
  for(size_t S=0;S<lim;S++) if(U c=dp[S]) { H ss=H(S), rest=full^ss; while(rest) { H bit=rest&(~rest+1); rest^=bit; int v=__builtin_ctzll(uint64_t(bit)); if(v==root || (g.adj[v]&ss)) { size_t T=S|size_t(bit); assert(UINT64_MAX-dp[T]>=c); dp[T]+=c; } } }
  return dp[lim-1];
}
int main(int ac,char**av) {
  if(ac==4 && std::string(av[1])=="exact") { int n=atoi(av[2]); if(n<2||n>4)return 2; G g(n); std::ofstream f(av[3]); f<<"n,root,x,y,h_root\n"; U sum=0; for(int r=0;r<g.N;r++) { U h=exact_root(g,r); sum+=h; f<<n<<','<<r<<','<<r%g.n<<','<<r/g.n<<','<<h<<'\n'; } f<<"# total_A1,"<<sum<<"\n"; return 0; }
  if(ac!=9) { std::cerr<<"usage: fixed_root_smc n root P R L threads out seed; or exact n out.csv\n"; return 2; }
  int n=atoi(av[1]),root=atoi(av[2]),P=atoi(av[3]),Rep=atoi(av[4]),L=atoi(av[5]),T=atoi(av[6]); std::string out=av[7]; U base=strtoull(av[8],nullptr,0);
  if(n<2||n>11||P<1||Rep<1||L<0||T<1||T>12||uint64_t(P)*uint64_t(n*n)>INT_MAX)return 2; G g(n); if(root<0||root>=g.N)return 2;
  std::ofstream meta(out+".metadata.json"); meta<<"{\"n\":"<<n<<",\"N\":"<<g.N<<",\"root\":"<<root<<",\"root_xy\":["<<root%g.n<<','<<root/g.n<<"],\"P\":"<<P<<",\"R\":"<<Rep<<",\"L\":"<<L<<",\"threads\":"<<T<<",\"seed_base\":\""<<base<<"\",\"normalizer\":\"h_root_count_of_onepeak_label_orders_not_divided_by_N_factorial\",\"mh\":\"symmetric_swap_replace_endpoint_rejection_root_fixed\"}\n";
  std::ofstream f(out); f<<"n,N,root,rep,L,seed,Z,min_parent_ess,dead,accept,proposals\n"; std::mutex mu; std::vector<std::thread>th;
  for(int t=0;t<T;t++) th.emplace_back([&,t]{ for(int rep=t;rep<Rep;rep+=T) { U sd=base+U(root)*1000000007+U(rep)*1009+U(L)*1000003; R z=run(g,root,P,L,sd); std::lock_guard<std::mutex>lk(mu); f<<n<<','<<g.N<<','<<root<<','<<rep<<','<<L<<','<<sd<<','<<std::setprecision(17)<<z.z<<','<<z.miness<<','<<z.dead<<','<<z.acc<<','<<z.props<<'\n'; f.flush(); }});
  for(auto&q:th)q.join();
}
