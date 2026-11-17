// Fixed-root one-peak SMC, N<=256.  Vertices are listed high label to low.
// probZ estimates h_root/(N-1)!; logZ is accumulated separately.
#include <array>
#include <bitset>
#include <algorithm>
#include <cassert>
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
using U=uint64_t; constexpr int MAXN=256;
struct W{std::array<unsigned char,MAXN>a{};};
struct G{int n,N;std::array<std::bitset<MAXN>,MAXN>adj{};explicit G(int q):n(q),N(q*q){for(int y=0;y<n;y++)for(int x=0;x<n;x++){int k=y*n+x;for(auto [dx,dy]:{std::pair{1,0},std::pair{-1,0},std::pair{0,1},std::pair{0,-1}}){int X=x+dx,Y=y+dy;if(X>=0&&X<n&&Y>=0&&Y<n)adj[k].set(Y*n+X);}}}};
struct R{double prob=1,logz=0,miness=0,acc=0,props=0;int dead=0;uint64_t terminal_hash=0;};
static int uni(std::mt19937_64&r,int m){assert(m>0);U cut=(U(0)-U(m))%U(m);for(;;){U x=r();if(x>=cut)return int(x%U(m));}}
static bool valid(const W&w,int r,const G&g){std::bitset<MAXN> seen;for(int t=0;t<r;t++){int v=w.a[t];if(seen.test(v))return false;if(t&&!(g.adj[v]&seen).any())return false;seen.set(v);}return true;}
static int allowed(const W&w,int r,const G&g,unsigned char*out){std::bitset<MAXN> seen,bd;for(int t=0;t<r;t++){int v=w.a[t];seen.set(v);bd|=g.adj[v];}if(r==0){for(int v=0;v<g.N;v++)out[v]=v;return g.N;}bd&=~seen;int k=0;for(int v=0;v<g.N;v++)if(bd.test(v))out[k++]=v;return k;}
static int peaks(const W&w,const G&g){std::array<int,MAXN>rank{};for(int i=0;i<g.N;i++)rank[w.a[i]]=i;int p=0;for(int v=0;v<g.N;v++){bool ok=true;for(int z=0;z<g.N;z++)if(g.adj[v].test(z)&&rank[z]<rank[v])ok=false;if(ok)p++;}return p;}
static W mh(W w,int r,int root,int L,const G&g,std::mt19937_64&rg,R&o){for(int q=0;q<L;q++){o.props++;int i=uni(rg,r),v=uni(rg,g.N),j=-1;for(int t=0;t<r;t++)if(w.a[t]==v){j=t;break;}W y=w;if(j>=0)std::swap(y.a[i],y.a[j]);else y.a[i]=v;if(y.a[0]==root&&valid(y,r,g)){w=y;o.acc++;}}return w;}
static R run(const G&g,int root,int P,int L,U seed){std::mt19937_64 rg(seed);std::vector<W>x(P);for(auto&w:x)w.a[0]=root;R o;unsigned char al[MAXN];for(int r=2;r<=g.N;r++){std::vector<int>a(P);double sw=0,s2=0;for(int i=0;i<P;i++){assert(x[i].a[0]==root&&valid(x[i],r-1,g));a[i]=allowed(x[i],r-1,g,al);sw+=a[i];s2+=double(a[i])*a[i];if(!a[i])o.dead++;}o.miness=(r==2?sw*sw/s2:std::min(o.miness,sw*sw/s2));double f=(sw/P)/double(g.N-r+1);o.prob*=f;o.logz+=std::log(f);if(!sw){o.prob=0;o.logz=-INFINITY;return o;}if(r==g.N){for(int i=0;i<P;i++)if(a[i]){int k=allowed(x[i],r-1,g,al);W y=x[i];y.a[r-1]=al[uni(rg,k)];assert(y.a[0]==root&&valid(y,r,g)&&peaks(y,g)==1);for(int t=0;t<g.N;t++)o.terminal_hash=o.terminal_hash*1315423911u+y.a[t];}return o;}std::vector<U>cdf(P);U total=0;for(int i=0;i<P;i++){total+=U(a[i]);cdf[i]=total;}std::vector<W>nx;nx.reserve(P);for(int q=0;q<P;q++){U u=U(uni(rg,int(total)));int j=int(std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin());if(j>=P)j=P-1;int k=allowed(x[j],r-1,g,al);W y=x[j];y.a[r-1]=al[uni(rg,k)];y=mh(y,r,root,L,g,rg,o);nx.push_back(y);}x.swap(nx);}return o;}
static U exact(const G&g,int root){assert(g.N<=16);size_t z=size_t(1)<<g.N;std::vector<U>d(z);d[0]=1;for(size_t S=0;S<z;S++)if(d[S])for(int v=0;v<g.N;v++)if(!(S>>v&1)&&(v==root||[&](){for(int q=0;q<g.N;q++)if((S>>q&1)&&g.adj[v].test(q))return true;return false;}())){assert(UINT64_MAX-d[S|(size_t(1)<<v)]>=d[S]);d[S|(size_t(1)<<v)]+=d[S];}return d[z-1];}
int main(int ac,char**av){if(ac==4&&std::string(av[1])=="exact"){int n=atoi(av[2]);G g(n);if(n<2||n>4)return 2;std::ofstream f(av[3]);f<<"n,root,x,y,h_root\n";U s=0;for(int r=0;r<g.N;r++){U h=exact(g,r);s+=h;f<<n<<','<<r<<','<<r%g.n<<','<<r/g.n<<','<<h<<'\n';}f<<"# total_A1,"<<s<<'\n';return 0;}if(ac!=9)return 2;int n=atoi(av[1]),root=atoi(av[2]),P=atoi(av[3]),Rep=atoi(av[4]),L=atoi(av[5]),T=atoi(av[6]);std::string out=av[7];U base=strtoull(av[8],0,0);if(n<2||n>16||P<1||Rep<1||L<0||T<1||T>12||uint64_t(P)*uint64_t(n*n)>INT_MAX)return 2;G g(n);if(root<0||root>=g.N)return 2;std::ofstream m(out+".metadata.json");m<<"{\"n\":"<<n<<",\"N\":"<<g.N<<",\"root\":"<<root<<",\"root_xy\":["<<root%g.n<<','<<root/g.n<<"],\"P\":"<<P<<",\"R\":"<<Rep<<",\"L\":"<<L<<",\"threads\":"<<T<<",\"seed_base\":\""<<base<<"\",\"normalizer\":\"probZ_estimates_h_root_over_(N-1)_factorial; logZ_is_log_probZ; root_is_first_highest_labelled_peak\"}\n";std::ofstream f(out);f<<"n,N,root,rep,L,seed,probZ,logZ,min_parent_ess,dead,accept,proposals,terminal_hash\n";std::mutex mu;std::vector<std::thread>th;for(int t=0;t<T;t++)th.emplace_back([&,t]{for(int q=t;q<Rep;q+=T){U sd=base+U(root)*1000000007+U(q)*1009+U(L)*1000003;R z=run(g,root,P,L,sd);std::lock_guard<std::mutex>lk(mu);f<<n<<','<<g.N<<','<<root<<','<<q<<','<<L<<','<<sd<<','<<std::setprecision(17)<<z.prob<<','<<z.logz<<','<<z.miness<<','<<z.dead<<','<<z.acc<<','<<z.props<<','<<z.terminal_hash<<'\n';f.flush();}});for(auto&q:th)q.join();}
