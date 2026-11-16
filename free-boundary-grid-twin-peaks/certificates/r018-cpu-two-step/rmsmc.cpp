// Ordered-prefix, fully-adapted resample-move SMC reference. n<=8, double diagnostics.
#include <array>
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
using U=uint64_t;
static std::ofstream DBGFILE("rmsmc-mh.csv");
static std::ostream* DBG=std::getenv("RMSMC_TRACE")?&DBGFILE:nullptr;
struct W{std::array<unsigned char,64>a{};};
struct R{double z=0,d=0,miness=0,acc=0,props=0; int dead=0;};
struct G{int n,N; std::array<U,64> adj{};G(int q):n(q),N(q*q){for(int y=0;y<n;y++)for(int x=0;x<n;x++){int k=y*n+x;for(auto [dx,dy]:{std::pair{1,0},{-1,0},{0,1},{0,-1}}){int X=x+dx,Y=y+dy;if(X>=0&&X<n&&Y>=0&&Y<n)adj[k]|=U(1)<<(Y*n+X);}}}};
static thread_local U LAST=0; static thread_local unsigned REJ=0;
static int uni(std::mt19937_64&r,int m){U cut=(U(0)-U(m))%U(m);REJ=0;do{LAST=r();if(LAST>=cut)return int(LAST%U(m));REJ++;}while(true);}
static bool valid(const W&w,int r,int tau,const G&g){U seen=0;int b=0;for(int t=0;t<r;t++){int v=w.a[t];if(!(g.adj[v]&seen)){b++;if(!(t==0||(t+1)==tau))return false;}seen|=U(1)<<v;}return b==(r<tau?1:2);}
static int allowed(const W&w,int r,int tau,const G&g,unsigned char*out){U seen=0,bd=0;for(int t=0;t<r;t++){int v=w.a[t];seen|=U(1)<<v;bd|=g.adj[v];}U m;if(r==0)m=(g.N==64?~U(0):(U(1)<<g.N)-1);else{bd&=~seen;m=(r+1==tau)?(((g.N==64?~U(0):(U(1)<<g.N)-1)^seen^bd)):bd;}int k=0;for(int v=0;v<g.N;v++)if(m&(U(1)<<v))out[k++]=v;return k;}
static int e2(const W&w,const G&g){U s=0;int p[2],b=0;for(int t=0;t<g.N;t++){int v=w.a[t];if(!(g.adj[v]&s))p[b++]=v;s|=U(1)<<v;}int dx=p[0]%g.n-p[1]%g.n,dy=p[0]/g.n-p[1]/g.n;return dx*dx+dy*dy;}
static W mh(W w,int r,int tau,int L,const G&g,std::mt19937_64&rg,R&z){for(int q=0;q<L;q++){z.props+=2;W old=w,y=w;int is[2],vs[2];for(int h=0;h<2;h++){is[h]=uni(rg,r);vs[h]=uni(rg,g.N);int j=-1;for(int t=0;t<r;t++)if(y.a[t]==vs[h])j=t;if(j>=0)std::swap(y.a[is[h]],y.a[j]);else y.a[is[h]]=vs[h];}bool ok=valid(y,r,tau,g);if(ok){w=y;z.acc++;}if(DBG){*DBG<<"B,";for(int t=0;t<r;t++)*DBG<<(t?":":"")<<int(old.a[t]);*DBG<<","<<tau<<","<<r<<","<<q<<","<<is[0]<<","<<vs[0]<<","<<is[1]<<","<<vs[1]<<","<<ok<<",";for(int t=0;t<r;t++)*DBG<<(t?":":"")<<int(w.a[t]);*DBG<<"\n";}}return w;}
static R run(const G&g,int tau,int P,int L,U seed,std::ofstream*trace=nullptr){std::mt19937_64 rg(seed);std::vector<W>x(P);R o;o.z=1;unsigned char al[64];for(int r=1;r<=g.N;r++){std::vector<int>a(P);double sw=0,s2=0;for(int i=0;i<P;i++){a[i]=allowed(x[i],r-1,tau,g,al);sw+=a[i];s2+=double(a[i])*a[i];if(!a[i])o.dead++;}o.miness=(r==1?sw*sw/s2:std::min(o.miness,sw*sw/s2));double mz=sw/P;o.z*=mz;if(!sw){o.z=o.d=0;return o;}if(r==g.N){double sf=0;for(int i=0;i<P;i++)if(a[i]){allowed(x[i],r-1,tau,g,al);W y=x[i];y.a[r-1]=al[uni(rg,a[i])];sf+=a[i]*e2(y,g);}o.d=(o.z/mz)*sf/P;return o;}std::vector<U>c(P);U s=0;for(int i=0;i<P;i++){s+=U(a[i]);c[i]=s;}std::vector<W>nx;nx.reserve(P);for(int q=0;q<P;q++){U u=U(uni(rg,int(s)));int j=int(std::upper_bound(c.begin(),c.end(),u)-c.begin());if(j>=P)j=P-1;allowed(x[j],r-1,tau,g,al);W y=x[j];y.a[r-1]=al[uni(rg,a[j])];y=mh(y,r,tau,L,g,rg,o);if(trace){*trace<<r;for(int t=0;t<r;t++)*trace<<","<<int(y.a[t]);*trace<<"\n";}nx.push_back(y);}x.swap(nx);}return o;}
int main(int ac,char**av){if(ac!=8){std::cerr<<"usage: rmsmc n P R L threads out seed\n";return 2;}int n=atoi(av[1]),P=atoi(av[2]),Rep=atoi(av[3]),L=atoi(av[4]),T=atoi(av[5]);std::string out=av[6];U seed=strtoull(av[7],0,0);if(n<2||n>8||P<1||Rep<1||T<1||L<0||uint64_t(P)*uint64_t(n*n)>INT_MAX)return 2;if(DBG&&T!=1){std::cerr<<"TRACE requires one worker\n";return 2;}G g(n);std::ofstream f(out),trace(out+".trace");f<<"n,tau,rep,L,seed,Z,D,Q,min_parent_ess,dead,accept,proposals\n";std::mutex mu;std::vector<std::thread> th;for(int t=0;t<T;t++)th.emplace_back([&,t]{for(int tau=2;tau<=g.N;tau++)for(int r=t;r<Rep;r+=T){U sd=seed+U(tau)*1000003+U(r)*1009+U(L)*1000000007;std::ofstream*tp=(n==3&&P==1&&Rep==1&&L==0&&tau==2&&r==0)?&trace:nullptr;R z=run(g,tau,P,L,sd,tp);std::lock_guard<std::mutex>lk(mu);f<<n<<','<<tau<<','<<r<<','<<L<<','<<sd<<','<<std::setprecision(17)<<z.z<<','<<z.d<<','<<(z.z?z.d/z.z/(g.N):NAN)<<','<<z.miness<<','<<z.dead<<','<<z.acc<<','<<z.props<<"\n";}});for(auto&q:th)q.join();}
