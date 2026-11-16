// Ordered-prefix, fully-adapted resample-move SMC reference. n<=11, double diagnostics.
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
using U=uint64_t;
using H=__uint128_t;
static std::ofstream DBGFILE("rmsmc-mh.csv");
static std::ostream* DBG=std::getenv("RMSMC_TRACE")?&DBGFILE:nullptr;
struct W{std::array<unsigned char,121>a{};};
struct R{double z=0,d=0,f=0,g=0,miness=0,acc=0,props=0; std::array<double,11> cols{}; int dead=0;};
struct G{int n,N; std::array<H,121> adj{};G(int q):n(q),N(q*q){for(int y=0;y<n;y++)for(int x=0;x<n;x++){int k=y*n+x;for(auto [dx,dy]:{std::pair{1,0},{-1,0},{0,1},{0,-1}}){int X=x+dx,Y=y+dy;if(X>=0&&X<n&&Y>=0&&Y<n)adj[k]|=H(1)<<(Y*n+X);}}}};
static thread_local U LAST=0; static thread_local unsigned REJ=0;
static int uni(std::mt19937_64&r,int m){U cut=(U(0)-U(m))%U(m);REJ=0;do{LAST=r();if(LAST>=cut)return int(LAST%U(m));REJ++;}while(true);}
static bool valid(const W&w,int r,int tau,const G&g){H seen=0;int b=0;for(int t=0;t<r;t++){int v=w.a[t];if(!(g.adj[v]&seen)){b++;if(!(t==0||(t+1)==tau))return false;}seen|=H(1)<<v;}return b==(r<tau?1:2);}
static int allowed(const W&w,int r,int tau,const G&g,unsigned char*out){
 H seen=0,bd=0;for(int t=0;t<r;t++){int v=w.a[t];seen|=H(1)<<v;bd|=g.adj[v];}
 H all=(H(1)<<g.N)-1,m;
 if(r==0)m=all;else{bd&=~seen;m=(r+1==tau)?(all^seen^bd):bd;}
 int k=0;for(int v=0;v<g.N;v++)if(m&(H(1)<<v))out[k++]=v;return k;
}
struct Moments { int d,f,g,p,q; };
static Moments moments(const W&w,const G&grid){
 H seen=0;int roots[2],b=0;std::array<int,121> rank{};
 for(int t=0;t<grid.N;t++){int v=w.a[t];assert(v<grid.N && !(seen&(H(1)<<v)));rank[v]=t;
  if(!(grid.adj[v]&seen)){assert(b<2);roots[b++]=v;}seen|=H(1)<<v;}
 assert(b==2);
 // Independent local strict-peak detector on rank labels, without bit masks.
 std::array<int,2> check{};int cb=0;
 for(int v=0;v<grid.N;v++){int x=v%grid.n,y=v/grid.n;bool peak=true;
  for(auto [X,Y]:{std::pair{x-1,y},{x+1,y},{x,y-1},{x,y+1}})
   if(X>=0&&X<grid.n&&Y>=0&&Y<grid.n&&rank[Y*grid.n+X]<rank[v])peak=false;
  if(peak){assert(cb<2);check[cb++]=v;}}
 assert(cb==2 && ((check[0]==roots[0]&&check[1]==roots[1])||(check[0]==roots[1]&&check[1]==roots[0])));
 int p=roots[0],q=roots[1],dx=p%grid.n-q%grid.n,dy=p/grid.n-q/grid.n;
 int px=2*(p%grid.n)-(grid.n-1),py=2*(p/grid.n)-(grid.n-1),qx=2*(q%grid.n)-(grid.n-1),qy=2*(q/grid.n)-(grid.n-1);
 Moments z{dx*dx+dy*dy,px*px+py*py+qx*qx+qy*qy,px*qx+py*qy,p,q};assert(z.f-2*z.g==4*z.d);return z;
}
static W mh(W w,int r,int tau,int L,const G&g,std::mt19937_64&rg,R&z){for(int q=0;q<L;q++){z.props++;int i=uni(rg,r);U ri=LAST;unsigned xi=REJ;int v=uni(rg,g.N);U rv=LAST;unsigned xv=REJ;int j=-1;for(int t=0;t<r;t++)if(w.a[t]==v)j=t;W y=w;for(int t=0;t<r;t++)if(DBG)*DBG<<(t?":":"M,")<<int(w.a[t]);if(j>=0)std::swap(y.a[i],y.a[j]);else y.a[i]=v;bool ok=valid(y,r,tau,g);if(ok){w=y;z.acc++;}if(DBG){*DBG<<","<<tau<<","<<r<<","<<q<<","<<ri<<","<<xi<<","<<rv<<","<<xv<<","<<i<<","<<v<<","<<ok<<",";for(int t=0;t<r;t++)*DBG<<(t?":":"")<<int(w.a[t]);*DBG<<"\n";}}return w;}
static R run(const G&g,int tau,int P,int L,U seed,std::ofstream*trace=nullptr){std::mt19937_64 rg(seed);std::vector<W>x(P);R o;o.z=1;unsigned char al[121];for(int r=1;r<=g.N;r++){std::vector<int>a(P);double sw=0,s2=0;for(int i=0;i<P;i++){a[i]=allowed(x[i],r-1,tau,g,al);sw+=a[i];s2+=double(a[i])*a[i];if(!a[i])o.dead++;}o.miness=(r==1?sw*sw/s2:std::min(o.miness,sw*sw/s2));double mz=sw/P;o.z*=mz;if(!sw){o.z=o.d=0;return o;}if(r==g.N){double sd=0,sf=0,sg=0;std::array<double,11>sc{};
 for(int i=0;i<P;i++)if(a[i]){allowed(x[i],r-1,tau,g,al);W y=x[i];y.a[r-1]=al[uni(rg,a[i])];auto mm=moments(y,g);
  sd+=a[i]*mm.d;sf+=a[i]*mm.f;sg+=a[i]*mm.g;sc[mm.p%g.n]+=a[i];sc[mm.q%g.n]+=a[i];}
 double factor=(o.z/mz)/P;o.d=(o.z/mz)*sd/P;o.f=factor*sf;o.g=factor*sg;for(int j=0;j<g.n;j++)o.cols[j]=factor*sc[j];return o;}
 std::vector<U>c(P);U s=0;for(int i=0;i<P;i++){s+=U(a[i]);c[i]=s;}std::vector<W>nx;nx.reserve(P);for(int q=0;q<P;q++){U u=U(uni(rg,int(s)));int j=int(std::upper_bound(c.begin(),c.end(),u)-c.begin());if(j>=P)j=P-1;allowed(x[j],r-1,tau,g,al);W y=x[j];y.a[r-1]=al[uni(rg,a[j])];y=mh(y,r,tau,L,g,rg,o);if(trace){*trace<<r;for(int t=0;t<r;t++)*trace<<","<<int(y.a[t]);*trace<<"\n";}nx.push_back(y);}x.swap(nx);}return o;}
int main(int ac,char**av){if(ac!=8){std::cerr<<"usage: rmsmc n P R L threads out seed\n";return 2;}int n=atoi(av[1]),P=atoi(av[2]),Rep=atoi(av[3]),L=atoi(av[4]),T=atoi(av[5]);std::string out=av[6];U seed=strtoull(av[7],0,0);if(n<2||n>11||P<1||Rep<1||T<1||L<0||uint64_t(P)*uint64_t(n*n)>INT_MAX)return 2;if(DBG&&T!=1){std::cerr<<"TRACE requires one worker\n";return 2;}std::ofstream metadata(out+".metadata.json");metadata<<"{\"n\":"<<n<<",\"P\":"<<P<<",\"R\":"<<Rep<<",\"L\":"<<L<<",\"threads\":"<<T<<",\"seed_base\":\""<<seed<<"\",\"F_definition\":\"count weight times sum of squared doubled centered root coordinates\",\"G_definition\":\"count weight times dot of doubled centered roots\"}\n";metadata.close();G g(n);std::ofstream f(out),trace(out+".trace");f<<"n,tau,rep,L,seed,Z,D,Q,min_parent_ess,dead,accept,proposals,F,G";for(int j=0;j<n;j++)f<<",COL"<<j;f<<"\n";std::mutex mu;std::vector<std::thread> th;for(int t=0;t<T;t++)th.emplace_back([&,t]{for(int tau=2;tau<=g.N;tau++)for(int r=t;r<Rep;r+=T){U sd=seed+U(tau)*1000003+U(r)*1009+U(L)*1000000007;std::ofstream*tp=(n==3&&P==1&&Rep==1&&L==0&&tau==2&&r==0)?&trace:nullptr;R z=run(g,tau,P,L,sd,tp);std::lock_guard<std::mutex>lk(mu);f<<n<<','<<tau<<','<<r<<','<<L<<','<<sd<<','<<std::setprecision(17)<<z.z<<','<<z.d<<','<<(z.z?z.d/z.z/(g.N):NAN)<<','<<z.miness<<','<<z.dead<<','<<z.acc<<','<<z.props<<','<<z.f<<','<<z.g;for(int j=0;j<n;j++)f<<','<<z.cols[j];f<<"\n";f.flush();}});for(auto&q:th)q.join();}
