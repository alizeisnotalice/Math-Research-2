// Exact first-merge root-depth ledger.  The depth is d_G(root,V\component).
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
using U=__uint128_t;
static std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static uint64_t choose64(int n,int k){if(k>n-k)k=n-k;uint64_t z=1;for(int i=1;i<=k;i++)z=z*(n-k+i)/i;return z;}
static uint32_t nextComb(uint32_t x){uint32_t c=x&-x,r=x+c;return (((r^x)>>2)/c)|r;}
static std::vector<uint32_t> masksOfSize(int n,int k){int N=n*n;std::vector<uint32_t>z;z.reserve(size_t(choose64(N,k)));if(k==N){z.push_back((uint32_t(1)<<N)-1);return z;}uint32_t x=(uint32_t(1)<<k)-1,cap=uint32_t(1)<<N;while(x<cap){z.push_back(x);x=nextComb(x);}return z;}
template<class F> static void parallelRanges(size_t c,int workers,F f){int used=std::min<int>(workers,int(c));std::vector<std::thread>ts;ts.reserve(used);for(int t=0;t<used;t++){size_t lo=c*size_t(t)/used,hi=c*size_t(t+1)/used;ts.emplace_back([=,&f](){f(lo,hi,t);});}for(auto&t:ts)t.join();}
static uint32_t component(uint32_t seed,uint32_t allowed,const std::array<uint32_t,25>&adj){uint32_t seen=seed,fr=seed;while(fr){uint32_t bit=fr&-fr;fr-=bit;uint32_t add=adj[__builtin_ctz(bit)]&allowed&~seen;seen|=add;fr|=add;}return seen;}
struct Orbit{int rep;std::vector<int> points;};
static std::vector<Orbit> d4orbits(int n){int N=n*n;std::vector<int>owner(N,-1);std::vector<Orbit>out;for(int v=0;v<N;v++)if(owner[v]<0){int x=v%n,y=v/n;std::array<std::array<int,2>,8> p{{{{x,y}},{{n-1-x,y}},{{x,n-1-y}},{{n-1-x,n-1-y}},{{y,x}},{{n-1-y,x}},{{y,n-1-x}},{{n-1-y,n-1-x}}}};std::vector<int>q;for(auto z:p){int w=z[1]*n+z[0];if(std::find(q.begin(),q.end(),w)==q.end())q.push_back(w);}int id=out.size();for(int w:q)owner[w]=id;out.push_back({*std::min_element(q.begin(),q.end()),q});}return out;}
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,workers=std::stoi(argv[3]);if(n<2||n>5||workers<1)return 2;
 uint32_t cap=uint32_t(1)<<N,full=cap-1;std::array<uint32_t,25>adj{};
 for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto alloc=[&](){return static_cast<U*>(std::calloc(cap,sizeof(U)));};U*a=alloc(),*h=alloc(),*ar=alloc();if(!a||!h||!ar){std::cerr<<"allocation failure\n";return 3;}
 auto start=std::chrono::steady_clock::now();
 for(int v=0;v<N;v++)a[uint32_t(1)<<v]=1;
 for(int k=2;k<=N;k++){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];U z=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int v=__builtin_ctz(bit);uint32_t T=S^bit;if(adj[v]&T)z+=a[T];}a[S]=z;}});}
 h[full]=1;for(int k=N-1;k>=1;k--){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i],out=full^S;U z=0;for(uint32_t q=out;q;q&=q-1){uint32_t bit=q&-q;if(adj[__builtin_ctz(bit)]&S)z+=h[S|bit];}h[S]=z;}});}
 const int MAXD=2*(n-1);std::vector<U>allHist(MAXD+1);auto orbits=d4orbits(n);
 std::ofstream o(argv[2]);o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers<<",\n  \"state_array_bytes\": "<<uint64_t(3)*cap*sizeof(U)<<",\n  \"depth_definition\": \"d_G(root,V\\\\component) in original grid\",\n  \"orbits\": [\n";
 bool ofirst=true;
 for(const auto&orb:orbits){int r=orb.rep;std::memset(ar,0,size_t(cap)*sizeof(U));ar[uint32_t(1)<<r]=1;
  for(int k=2;k<=N;k++){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];if(!(S&(uint32_t(1)<<r)))continue;U z=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int v=__builtin_ctz(bit);uint32_t T=S^bit;if(adj[v]&T)z+=ar[T];}ar[S]=z;}});}
  std::array<uint32_t,9> ball{};for(int d=0;d<=MAXD;d++)for(int v=0;v<N;v++){int x=v%n,y=v/n,rx=r%n,ry=r/n;if(std::abs(x-rx)+std::abs(y-ry)<=d)ball[d]|=uint32_t(1)<<v;}
  std::vector<std::vector<U>> local{size_t(workers),std::vector<U>(MAXD+1)};
  for(int s=3;s<=N;s++){auto masks=masksOfSize(n,s);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int tid){auto&hist=local[tid];for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];if(!h[S]||!a[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int w=__builtin_ctz(bit);uint32_t T=S^bit,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;U base=h[S]*U(choose64(s-1,__builtin_popcount(A)));if(A&(uint32_t(1)<<r)){int d=1;while(d<=MAXD&&(ball[d-1]&~A)==0)d++;hist[d-1]+=base*ar[A]*a[B];}if(B&(uint32_t(1)<<r)){int d=1;while(d<=MAXD&&(ball[d-1]&~B)==0)d++;hist[d-1]+=base*a[A]*ar[B];}}}});}
  std::vector<U> one(MAXD+1);for(auto&v:local)for(int d=1;d<=MAXD;d++)one[d]+=v[d];U oneTotal=0;for(int d=1;d<=MAXD;d++){oneTotal+=one[d];allHist[d]+=U(orb.points.size())*one[d];}
  if(!ofirst)o<<",\n";ofirst=false;o<<"    {\"representative\": "<<r<<", \"orbit_size\": "<<orb.points.size()<<", \"one_root_total\": \""<<dec(oneTotal)<<"\", \"histogram\": {";bool hf=true;for(int d=1;d<=MAXD;d++)if(one[d]){if(!hf)o<<", ";hf=false;o<<"\""<<d<<"\": \""<<dec(one[d])<<"\"";}o<<"}}";
 }
 U allTotal=0,second=0;for(int d=1;d<=MAXD;d++){allTotal+=allHist[d];second+=U(d)*d*allHist[d];}
 o<<"\n  ],\n  \"all_roots_histogram\": {";bool hf=true;for(int d=1;d<=MAXD;d++)if(allHist[d]){if(!hf)o<<", ";hf=false;o<<"\""<<d<<"\": \""<<dec(allHist[d])<<"\"";}o<<"},\n  \"all_roots_total\": \""<<dec(allTotal)<<"\",\n  \"all_roots_depth_squared_sum\": \""<<dec(second)<<"\",\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";
 o.close();std::cout<<"n="<<n<<" root_total="<<dec(allTotal)<<" depth2="<<dec(second)<<" orbits="<<orbits.size()<<"\n";std::free(a);std::free(h);std::free(ar);
}
