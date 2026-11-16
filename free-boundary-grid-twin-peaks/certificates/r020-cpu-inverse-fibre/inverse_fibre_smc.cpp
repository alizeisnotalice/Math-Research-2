// Ordered-prefix fully-adapted resample-move SMC, restricted to ONE peak.
// Proposal, integer CDF resampling, RNG rejection sampling and MH are copied
// unchanged from r019-cpu-root-moments/covsmc.cpp.  tau=N+1 means no second
// birth.  Terminal words are independently checked for exactly one peak.
#include <array>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <thread>
#include <mutex>
#include <algorithm>
#include <sstream>
using U=uint64_t; using H=__uint128_t;
static std::ofstream DBGFILE("rmsmc-mh.csv");
static std::ostream* DBG=std::getenv("RMSMC_TRACE")?&DBGFILE:nullptr;
struct W{std::array<unsigned char,121>a{};};
struct R{double z=0,m2z=0,m3z=0,miness=0,acc=0,props=0;int dead=0;};
struct G{int n,N;std::array<H,121>adj{};G(int q):n(q),N(q*q){for(int y=0;y<n;y++)for(int x=0;x<n;x++){int k=y*n+x;for(auto[dx,dy]:{std::pair{1,0},{-1,0},{0,1},{0,-1}}){int X=x+dx,Y=y+dy;if(X>=0&&X<n&&Y>=0&&Y<n)adj[k]|=H(1)<<(Y*n+X);}}}};
static thread_local U LAST=0;static thread_local unsigned REJ=0;
static int uni(std::mt19937_64&r,int m){U cut=(U(0)-U(m))%U(m);REJ=0;do{LAST=r();if(LAST>=cut)return int(LAST%U(m));REJ++;}while(true);}
static bool valid(const W&w,int r,int tau,const G&g){H seen=0;int b=0;for(int t=0;t<r;t++){int v=w.a[t];if(!(g.adj[v]&seen)){b++;if(!(t==0||(t+1)==tau))return false;}seen|=H(1)<<v;}return b==(r<tau?1:2);}
static int allowed(const W&w,int r,int tau,const G&g,unsigned char*out){H seen=0,bd=0;for(int t=0;t<r;t++){int v=w.a[t];seen|=H(1)<<v;bd|=g.adj[v];}H all=(H(1)<<g.N)-1,m;if(r==0)m=all;else{bd&=~seen;m=(r+1==tau)?(all^seen^bd):bd;}int k=0;for(int v=0;v<g.N;v++)if(m&(H(1)<<v))out[k++]=v;return k;}
static int md(int a,int b,const G&g){return std::abs(a%g.n-b%g.n)+std::abs(a/g.n-b/g.n);}
static std::vector<int> canonical_path(int u,int v,const G&g){std::vector<int>p{u};while(u!=v){int rem=md(u,v,g),nx=g.N;for(int x=0;x<g.N;x++)if((g.adj[u]&(H(1)<<x))&&md(x,v,g)==rem-1)nx=std::min(nx,x);assert(nx<g.N);p.push_back(nx);u=nx;}return p;}
static int independent_one_peak_check(const W&w,const G&g){std::array<int,121>rank{};for(int t=0;t<g.N;t++)rank[w.a[t]]=t;int b=0;for(int v=0;v<g.N;v++){bool peak=true;for(int z=0;z<g.N;z++)if((g.adj[v]&(H(1)<<z))&&rank[z]<rank[v])peak=false;if(peak)b++;}return b;}
static U checked_add(U a,U b){assert(UINT64_MAX-a>=b);return a+b;}
// Exact fibre size of the canonical greedy shortest-path lift, for d=2 or 3.
static U inverse_fibre(const W&sigma,int d,const G&g){
  assert(independent_one_peak_check(sigma,g)==1);int u=sigma.a[0],plen=d-1;std::array<int,2>P{};for(int j=0;j<plen;j++)P[j]=sigma.a[1+j];
  std::vector<int>pref{u};for(int j=0;j<plen;j++)pref.push_back(P[j]);std::vector<int>tau{u};for(int j=d;j<g.N;j++)tau.push_back(sigma.a[j]);
  H pm=0;for(int x:pref)pm|=H(1)<<x;std::vector<H>prefix(tau.size()+1);for(size_t i=0;i<tau.size();i++)prefix[i+1]=prefix[i]|(H(1)<<tau[i]);
  std::array<H,4> inserted{};for(int s=1;s<(1<<plen);s++){int bit=__builtin_ctz(unsigned(s));inserted[s]=inserted[s^(1<<bit)]|(H(1)<<P[bit]);}
  U answer=0;for(int v=0;v<g.N;v++){
    if((pm&(H(1)<<v))||md(u,v,g)!=d||!(g.adj[P[plen-1]]&(H(1)<<v)))continue;
    auto path=canonical_path(u,v,g);bool same=int(path.size())==d+1;for(int j=0;j<d&&same;j++)same=path[j]==pref[j];if(!same)continue;
    const int states=(int(tau.size())+1)*(1<<plen),full=(1<<plen)-1;std::vector<U>dp(states),nd(states);auto at=[&](int i,int s){return i*(1<<plen)+s;};dp[at(1,0)]=1;
    for(int step=1;step<g.N;step++){std::fill(nd.begin(),nd.end(),0);for(int i=0;i<=int(tau.size());i++)for(int s=0;s<=full;s++){U c=dp[at(i,s)];if(!c)continue;H placed=prefix[i]|inserted[s];for(int j=0;j<plen;j++)if(!(s&(1<<j))&&(g.adj[P[j]]&placed))nd[at(i,s|(1<<j))]=checked_add(nd[at(i,s|(1<<j))],c);if(i<int(tau.size())){int x=tau[i];bool born=!(g.adj[x]&placed);if((x==v&&born)||(x!=v&&!born))nd[at(i+1,s)]=checked_add(nd[at(i+1,s)],c);}}
      dp.swap(nd);
    }
    answer=checked_add(answer,dp[at(int(tau.size()),full)]);
  }
  return answer;
}
static std::pair<U,U> fibres(const W&w,const G&g){return {inverse_fibre(w,2,g),inverse_fibre(w,3,g)};}
// Unchanged ordered-word MH proposal/rejection kernel.
static W mh(W w,int r,int tau,int L,const G&g,std::mt19937_64&rg,R&z){for(int q=0;q<L;q++){z.props++;int i=uni(rg,r);U ri=LAST;unsigned xi=REJ;int v=uni(rg,g.N);U rv=LAST;unsigned xv=REJ;int j=-1;for(int t=0;t<r;t++)if(w.a[t]==v)j=t;W y=w;for(int t=0;t<r;t++)if(DBG)*DBG<<(t?":":"M,")<<int(w.a[t]);if(j>=0)std::swap(y.a[i],y.a[j]);else y.a[i]=v;bool ok=valid(y,r,tau,g);if(ok){w=y;z.acc++;}if(DBG){*DBG<<","<<tau<<","<<r<<","<<q<<","<<ri<<","<<xi<<","<<rv<<","<<xv<<","<<i<<","<<v<<","<<ok<<",";for(int t=0;t<r;t++)*DBG<<(t?":":"")<<int(w.a[t]);*DBG<<"\n";}}return w;}
static R run(const G&g,int P,int L,U seed,std::ofstream*trace=nullptr){const int tau=g.N+1;std::mt19937_64 rg(seed);std::vector<W>x(P);R o;o.z=1;unsigned char al[121];for(int r=1;r<=g.N;r++){std::vector<int>a(P);double sw=0,s2=0;for(int i=0;i<P;i++){a[i]=allowed(x[i],r-1,tau,g,al);sw+=a[i];s2+=double(a[i])*a[i];if(!a[i])o.dead++;}o.miness=(r==1?sw*sw/s2:std::min(o.miness,sw*sw/s2));double mz=sw/P;o.z*=mz;if(!sw){o.z=o.m2z=o.m3z=0;return o;}if(r==g.N){double sm2=0,sm3=0;for(int i=0;i<P;i++)if(a[i]){allowed(x[i],r-1,tau,g,al);W y=x[i];y.a[r-1]=al[uni(rg,a[i])];assert(independent_one_peak_check(y,g)==1);auto[m2,m3]=fibres(y,g);sm2+=a[i]*double(m2);sm3+=a[i]*double(m3);}double factor=(o.z/mz)/P;o.m2z=factor*sm2;o.m3z=factor*sm3;return o;}
 std::vector<U>c(P);U s=0;for(int i=0;i<P;i++){s+=U(a[i]);c[i]=s;}std::vector<W>nx;nx.reserve(P);for(int q=0;q<P;q++){U u=U(uni(rg,int(s)));int j=int(std::upper_bound(c.begin(),c.end(),u)-c.begin());if(j>=P)j=P-1;allowed(x[j],r-1,tau,g,al);W y=x[j];y.a[r-1]=al[uni(rg,a[j])];y=mh(y,r,tau,L,g,rg,o);if(trace){*trace<<r;for(int t=0;t<r;t++)*trace<<","<<int(y.a[t]);*trace<<"\n";}nx.push_back(y);}x.swap(nx);}return o;}
static std::string key(const W&w,const G&g){std::string s;for(int i=0;i<g.N;i++){if(i)s+='-';s+=std::to_string(int(w.a[i]));}return s;}
static void audit_rec(W&w,int r,const G&g,std::ofstream&out,U&count){unsigned char al[121];int k=allowed(w,r,g.N+1,g,al);for(int j=0;j<k;j++){w.a[r]=al[j];if(r+1==g.N){assert(independent_one_peak_check(w,g)==1);auto[m2,m3]=fibres(w,g);out<<key(w,g)<<","<<m2<<","<<m3<<"\n";count++;}else audit_rec(w,r+1,g,out,count);}}
int main(int ac,char**av){
  if(ac==3&&std::string(av[1])=="audit-n3"){G g(3);std::ofstream o(av[2]);o<<"word,m2,m3\n";W w;U count=0;audit_rec(w,0,g,o,count);assert(count==29568);return 0;}
  if(ac!=8){std::cerr<<"usage: inverse_fibre_smc n P R L threads out seed; or audit-n3 out.csv\n";return 2;}
  int n=atoi(av[1]),P=atoi(av[2]),Rep=atoi(av[3]),L=atoi(av[4]),T=atoi(av[5]);std::string out=av[6];U seed=strtoull(av[7],0,0);
  if(n<2||n>11||P<1||Rep<1||T<1||T>8||L<0||uint64_t(P)*uint64_t(n*n)>INT_MAX)return 2;if(DBG&&T!=1){std::cerr<<"TRACE requires one worker\n";return 2;}
  G g(n);const int tau=g.N+1;std::ofstream metadata(out+".metadata.json");metadata<<"{\"n\":"<<n<<",\"P\":"<<P<<",\"R\":"<<Rep<<",\"L\":"<<L<<",\"threads\":"<<T<<",\"seed_base\":\""<<seed<<"\",\"tau\":"<<tau<<",\"target\":\"A2_distance_over_A1_equals_E_onepeak_m_distance\",\"canonical_rule\":\"least-numbered neighbor reducing Manhattan distance\"}\n";metadata.close();
  std::ofstream f(out),trace(out+".trace");f<<"n,tau,rep,L,seed,Z,Zmean_m2,Zmean_m3,weighted_ratio_m2,weighted_ratio_m3,min_parent_ess,dead,accept,proposals\n";std::mutex mu;std::vector<std::thread>th;for(int t=0;t<T;t++)th.emplace_back([&,t]{for(int r=t;r<Rep;r+=T){U sd=seed+U(tau)*1000003+U(r)*1009+U(L)*1000000007;std::ofstream*tp=(n==3&&P==1&&Rep==1&&L==0&&r==0)?&trace:nullptr;R z=run(g,P,L,sd,tp);std::lock_guard<std::mutex>lk(mu);f<<n<<','<<tau<<','<<r<<','<<L<<','<<sd<<','<<std::setprecision(17)<<z.z<<','<<z.m2z<<','<<z.m3z<<','<<(z.z?z.m2z/z.z:NAN)<<','<<(z.z?z.m3z/z.z:NAN)<<','<<z.miness<<','<<z.dead<<','<<z.acc<<','<<z.props<<"\n";f.flush();}});for(auto&q:th)q.join();
}
