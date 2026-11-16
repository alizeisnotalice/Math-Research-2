// Multithreaded CPU sequential importance sampler for the original K=2 grid law.
// Numerical diagnostic only: ideal exact-uniform sampling gives unbiased count
// numerators; floating output and their ratios are not probability certificates.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

constexpr int MW=4; // up to 256 vertices
using U64=uint64_t;
struct Bits { std::array<U64,MW> w{}; };
static Bits operator|(Bits a,const Bits&b){for(int i=0;i<MW;i++)a.w[i]|=b.w[i];return a;}
static Bits operator&(Bits a,const Bits&b){for(int i=0;i<MW;i++)a.w[i]&=b.w[i];return a;}
static Bits diff(Bits a,const Bits&b){for(int i=0;i<MW;i++)a.w[i]&=~b.w[i];return a;}
static int pop(const Bits&a){int z=0;for(U64 x:a.w)z+=__builtin_popcountll(x);return z;}
static bool has(const Bits&a,int v){return (a.w[v>>6]>>(v&63))&1;}
static void set(Bits&a,int v){a.w[v>>6]|=U64(1)<<(v&63);}
struct RNG { std::mt19937_64 e; explicit RNG(U64 s):e(s){}; int below(int n){ U64 lim=(std::numeric_limits<U64>::max()/U64(n))*U64(n); U64 x; do{x=e();}while(x>=lim); return int(x%U64(n)); }};
static int kth(const Bits&a,int k){for(int q=0;q<MW;q++){U64 x=a.w[q];int c=__builtin_popcountll(x);if(k>=c){k-=c;continue;}while(k--){x&=x-1;}return q*64+__builtin_ctzll(x);}throw std::runtime_error("bad kth");}
struct Grid { int n,N; Bits full; std::vector<Bits> adj; Grid(int nn):n(nn),N(nn*nn),adj(N){if(n<2||N>256)throw std::runtime_error("n must be 2..16");for(int v=0;v<N;v++)set(full,v);for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:{std::pair{1,0},std::pair{-1,0},std::pair{0,1},std::pair{0,-1}}){int X=x+d.first,Y=y+d.second;if(X>=0&&X<n&&Y>=0&&Y<n)set(adj[v],Y*n+X);}}} };
struct Stat { long double a1=0,a1sq=0,a2=0,a2sq=0,e2=0,rho2=0,farhalf=0,farquarter=0,maxa1=0,maxa2=0; U64 m=0,k2=0,dead=0; };
static void add(Stat&a,const Stat&b){a.a1+=b.a1;a.a1sq+=b.a1sq;a.a2+=b.a2;a.a2sq+=b.a2sq;a.e2+=b.e2;a.rho2+=b.rho2;a.farhalf+=b.farhalf;a.farquarter+=b.farquarter;a.maxa1=std::max(a.maxa1,b.maxa1);a.maxa2=std::max(a.maxa2,b.maxa2);a.m+=b.m;a.k2+=b.k2;a.dead+=b.dead;}
struct Sample { std::vector<int> order; std::array<int,2> peaks; int rho; };
static std::string num(long double x){std::ostringstream s;s<<std::setprecision(21)<<std::scientific<<x;return s.str();}
static void oneK2(const Grid&g,RNG&r,Stat&z,std::vector<Sample>& keep){
 Bits S{},front{}; std::array<int,2> peaks{-1,-1}; std::vector<int> ord;ord.reserve(g.N);int born=0;long double W=1;
 for(int t=0;t<g.N;t++){
  Bits avail=born<2?diff(g.full,S):diff(front,S);int a=pop(avail);if(a<=0)throw std::runtime_error("empty allowed");W*=a;int v=kth(avail,r.below(a)); bool isborn=!pop(g.adj[v]&S);if(isborn){if(born>=2)throw std::runtime_error("third peak");peaks[born++]=v;}set(S,v);front=front|g.adj[v];ord.push_back(v);
 }
 z.m++;if(born!=2) return;z.k2++;z.a2+=W;z.a2sq+=W*W;z.maxa2=std::max(z.maxa2,W);int x1=peaks[0]%g.n,y1=peaks[0]/g.n,x2=peaks[1]%g.n,y2=peaks[1]/g.n;int dx=std::abs(x1-x2),dy=std::abs(y1-y2),e=dx*dx+dy*dy,rho=dx+dy;z.e2+=W*e;z.rho2+=W*rho*rho;if(rho>= (g.n+1)/2)z.farhalf+=W;if(4*rho>g.n)z.farquarter+=W;if(keep.size()<5)keep.push_back({ord,peaks,rho});
}
static void oneK2stratum(const Grid&g,int tau,RNG&r,Stat&z,std::vector<Sample>&keep){
 Bits S{},front{};std::array<int,2>peaks{-1,-1};std::vector<int>ord;ord.reserve(g.N);long double W=1;z.m++;
 for(int t=0;t<g.N;t++){
  int rank=t+1;Bits un=diff(g.full,S),avail;
  if(rank==1)avail=un;else if(rank<tau)avail=diff(front,S);else if(rank==tau)avail=diff(un,diff(front,S));else avail=diff(front,S);
  int a=pop(avail);if(!a){z.dead++;return;}W*=a;int v=kth(avail,r.below(a));bool born=!pop(g.adj[v]&S);
  if(rank==1){if(!born)throw std::runtime_error("first must birth");peaks[0]=v;}
  else if(rank<tau){if(born)throw std::runtime_error("forced boundary birthed");}
  else if(rank==tau){if(!born)throw std::runtime_error("forced nonboundary did not birth");peaks[1]=v;}
  else if(born)throw std::runtime_error("post-tau third birth");
  set(S,v);front=front|g.adj[v];ord.push_back(v);
 }
 z.k2++;z.a2+=W;z.a2sq+=W*W;z.maxa2=std::max(z.maxa2,W);int x1=peaks[0]%g.n,y1=peaks[0]/g.n,x2=peaks[1]%g.n,y2=peaks[1]/g.n,dx=std::abs(x1-x2),dy=std::abs(y1-y2),e=dx*dx+dy*dy,rho=dx+dy;z.e2+=W*e;z.rho2+=W*rho*rho;if(rho>=(g.n+1)/2)z.farhalf+=W;if(4*rho>g.n)z.farquarter+=W;if(keep.size()<5)keep.push_back({ord,peaks,rho});
}
static void oneK1(const Grid&g,RNG&r,Stat&z){Bits S{},front{};long double W=1;for(int t=0;t<g.N;t++){Bits avail=t==0?diff(g.full,S):diff(front,S);int a=pop(avail);if(a<=0)throw std::runtime_error("empty K1 allowed");W*=a;int v=kth(avail,r.below(a));set(S,v);front=front|g.adj[v];}z.a1+=W;z.a1sq+=W*W;z.maxa1=std::max(z.maxa1,W);}
static long double ess(long double s,long double s2){return s2? s*s/s2:0;}
static long double ratio(long double x,long double y){return y?x/y:0;}
static void write_batch(std::ofstream&o,const Grid&g,int tau,int b,const Stat&z,double sec,U64 seed){
 long double M=z.m;
 o<<"{\"kind\":\"batch\",\"n\":"<<g.n<<",\"tau\":"<<tau<<",\"batch\":"<<b
  <<",\"samples\":"<<z.m<<",\"raw_K2\":"<<z.k2<<",\"dead\":"<<z.dead<<",\"seed\":"<<seed
  <<",\"seconds\":"<<std::setprecision(10)<<sec
  <<",\"A1_mean\":\""<<num(z.a1/M)<<"\",\"A1_ess\":\""<<num(ess(z.a1,z.a1sq))
  <<"\",\"A1_top_share\":\""<<num(ratio(z.maxa1,z.a1))
  <<"\",\"A2_mean\":\""<<num(z.a2/M)<<"\",\"A2_ess\":\""<<num(ess(z.a2,z.a2sq))
  <<"\",\"A2_top_share\":\""<<num(ratio(z.maxa2,z.a2))
  <<"\",\"Q_e2\":\""<<num(ratio(z.e2,g.N*z.a2))
  <<"\",\"rho2_mean\":\""<<num(ratio(z.rho2,z.a2))
  <<"\",\"far_half_weight\":\""<<num(ratio(z.farhalf,z.a2))
  <<"\",\"far_quarter_weight\":\""<<num(ratio(z.farquarter,z.a2))<<"\"}"<<'\n';
}
static U64 splitmix(U64 x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
static void write_samples(const std::string& path,int n,const std::vector<Sample>& ss){
 std::ofstream o(path);o<<"{\"n\":"<<n<<",\"samples\":[";
 for(size_t i=0;i<ss.size();i++){if(i)o<<',';o<<"{\"descending_order\":[";for(size_t j=0;j<ss[i].order.size();j++){if(j)o<<',';o<<ss[i].order[j];}o<<"],\"reported_peaks\":["<<ss[i].peaks[0]<<','<<ss[i].peaks[1]<<"],\"reported_rho\":"<<ss[i].rho<<'}';}o<<"]}\n";
}
int main(int argc,char**argv){try{
 if(argc!=6){std::cerr<<"usage: original_law_sis n batches samples_per_tau_batch threads output_dir\n";return 2;}
 int n=std::stoi(argv[1]), B=std::stoi(argv[2]), per=std::stoi(argv[3]), nt=std::stoi(argv[4]);std::string out=argv[5];
 if(B<=0||per<=0||nt<=0)throw std::runtime_error("positive batches/samples/threads required");
 if(n>8)throw std::runtime_error("direct W and W^2 are restricted to n<=8; use a log-sum-exp implementation for n>8");
 Grid g(n);std::string cmd="mkdir -p '"+out+"'";if(std::system(cmd.c_str()))throw std::runtime_error("cannot create output");
 std::ofstream js(out+"/batches.jsonl");js<<"{\"kind\":\"metadata\",\"n\":"<<n<<",\"batches\":"<<B<<",\"samples_per_tau_batch\":"<<per<<",\"threads\":"<<nt<<",\"generator\":\"mt19937_64 with rejection-uniform index\",\"proposal_K2\":\"forced second birth rank tau\",\"proposal_K1\":\"first:any; then frontier\",\"numerical_only\":true}\n";
 std::ofstream su(out+"/summary.jsonl");std::vector<Sample> saved;
 for(int tau=2;tau<=g.N;tau++){Stat all;
 for(int b=0;b<B;b++){
  auto start=std::chrono::steady_clock::now();std::vector<Stat> local(nt);std::vector<std::vector<Sample>> keep(nt);std::vector<std::thread> th;
  U64 seed=splitmix(0xA022C001ULL+U64(n)*1000003ULL+U64(tau)*1009ULL+U64(b));
  for(int t=0;t<nt;t++)th.emplace_back([&,t]{RNG r(splitmix(seed+U64(t)));int lo=per*t/nt,hi=per*(t+1)/nt;for(int i=lo;i<hi;i++){oneK2stratum(g,tau,r,local[t],keep[t]);if(tau==2)oneK1(g,r,local[t]);}});
  for(auto&x:th)x.join();Stat z;for(int t=0;t<nt;t++){add(z,local[t]);for(auto&s:keep[t])if(saved.size()<12)saved.push_back(std::move(s));}add(all,z);double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();write_batch(js,g,tau,b,z,sec,seed);js.flush();std::cerr<<"n="<<n<<" tau="<<tau<<" batch="<<b<<" samples="<<z.m<<" K2="<<z.k2<<" dead="<<z.dead<<" sec="<<sec<<" A2ESS="<<(double)ess(z.a2,z.a2sq)<<"\n";
 }
 long double M=all.m;su<<"{\"kind\":\"summary\",\"n\":"<<n<<",\"tau\":"<<tau<<",\"samples\":"<<all.m<<",\"dead\":"<<all.dead<<",\"A1_mean\":\""<<num(all.a1/M)<<"\",\"A1_ess\":\""<<num(ess(all.a1,all.a1sq))<<"\",\"A2_mean\":\""<<num(all.a2/M)<<"\",\"A2_ess\":\""<<num(ess(all.a2,all.a2sq))<<"\",\"A2_top_share\":\""<<num(ratio(all.maxa2,all.a2))<<"\",\"Q_e2\":\""<<num(ratio(all.e2,g.N*all.a2))<<"\",\"rho2_mean\":\""<<num(ratio(all.rho2,all.a2))<<"\",\"far_half_weight\":\""<<num(ratio(all.farhalf,all.a2))<<"\",\"far_quarter_weight\":\""<<num(ratio(all.farquarter,all.a2))<<"\"}\n";
 }
 write_samples(out+"/saved-k2-labels.json",n,saved);
 return 0;
}catch(const std::exception&e){std::cerr<<"original_law_sis: "<<e.what()<<'\n';return 1;}}
