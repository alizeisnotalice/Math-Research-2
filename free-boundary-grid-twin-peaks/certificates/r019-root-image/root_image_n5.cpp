// Exact n<=5 root-image ledger for the canonical detachable-branch insertion.
// It computes sum W q_B by one specified-root subset DP per grid vertex.
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
static std::string dec(U x){if(!x)return "0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static uint64_t choose64(int n,int k){if(k>n-k)k=n-k;uint64_t z=1;for(int i=1;i<=k;i++)z=z*(n-k+i)/i;return z;}
template<class F> static void ranges(size_t count,int workers,F f){int used=std::min<int>(workers,int(count));std::vector<std::thread>ts;ts.reserve(used);for(int t=0;t<used;t++){size_t lo=count*size_t(t)/used,hi=count*size_t(t+1)/used;ts.emplace_back([=,&f]{f(lo,hi,t);});}for(auto&t:ts)t.join();}
static uint32_t component(uint32_t seed,uint32_t allowed,const std::array<uint32_t,25>& adj){uint32_t seen=seed,front=seed;while(front){uint32_t bit=front&-front;front-=bit;uint32_t add=adj[__builtin_ctz(bit)]&allowed&~seen;seen|=add;front|=add;}return seen;}
struct Agg{U w=0,img=0;uint64_t states=0;};
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,workers=std::stoi(argv[3]);if(n<3||n>5||workers<1||workers>4)return 2;
 uint32_t cap=uint32_t(1)<<N,full=cap-1;std::array<uint32_t,25>adj{};
 for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto alloc=[&](){return static_cast<U*>(std::calloc(cap,sizeof(U)));};U*a=alloc(),*h=alloc(),*ar=alloc();if(!a||!h||!ar)return 3;
 auto start=std::chrono::steady_clock::now();std::vector<std::vector<uint32_t>> layer(N+1);
 for(uint32_t S=1;S<cap;S++)layer[__builtin_popcount(S)].push_back(S);
 for(int v=0;v<N;v++)a[uint32_t(1)<<v]=1;
 for(int k=2;k<=N;k++)ranges(layer[k].size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=layer[k][i];U z=0;for(uint32_t q=S;q;q&=q-1){uint32_t b=q&-q,T=S^b;if(adj[__builtin_ctz(b)]&T)z+=a[T];}a[S]=z;}});
 h[full]=1;for(int k=N-1;k>=1;k--)ranges(layer[k].size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=layer[k][i],out=full^S;U z=0;for(uint32_t q=out;q;q&=q-1){uint32_t b=q&-q;if(adj[__builtin_ctz(b)]&S)z+=h[S|b];}h[S]=z;}});
 std::vector<Agg> byB(N+1),byRank(N+1);std::vector<Agg> allLocal{size_t(workers)};
 ranges(cap,workers,[&](size_t lo,size_t hi,int tid){Agg&z=allLocal[tid];for(uint32_t S=uint32_t(lo);S<uint32_t(hi);S++){int rank=__builtin_popcount(S);if(rank<3||!h[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t wb=q&-q,T=S^wb,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;int w=__builtin_ctz(wb);if(!(adj[w]&A)||!(adj[w]&B))continue;if(__builtin_ctz(A)>__builtin_ctz(B))std::swap(A,B);int aa=__builtin_popcount(A),bb=__builtin_popcount(B);U W=h[S]*U(choose64(aa+bb,aa))*a[A]*a[B];z.w+=W;z.states++;}}});
 Agg all{};for(auto z:allLocal){all.w+=z.w;all.states+=z.states;}
 // Re-scan states once to form W strata.  (No root counts are used here.)
 for(uint32_t S=1;S<cap;S++){int rank=__builtin_popcount(S);if(rank<3||!h[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t wb=q&-q,T=S^wb,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;int w=__builtin_ctz(wb);if(!(adj[w]&A)||!(adj[w]&B))continue;if(__builtin_ctz(A)>__builtin_ctz(B))std::swap(A,B);int aa=__builtin_popcount(A),bb=__builtin_popcount(B);U W=h[S]*U(choose64(aa+bb,aa))*a[A]*a[B];byB[bb].w+=W;byB[bb].states++;byRank[rank].w+=W;byRank[rank].states++;}}
 // Root r is summed only on the canonical B side and only when r~w.
 for(int r=0;r<N;r++){
  std::memset(ar,0,size_t(cap)*sizeof(U));ar[uint32_t(1)<<r]=1;
  for(int k=2;k<=N;k++)ranges(layer[k].size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=layer[k][i];if(!(S&(uint32_t(1)<<r)))continue;U z=0;for(uint32_t q=S;q;q&=q-1){uint32_t b=q&-q,T=S^b;if(adj[__builtin_ctz(b)]&T)z+=ar[T];}ar[S]=z;}});
  std::vector<std::vector<Agg>> localB{size_t(workers),std::vector<Agg>(N+1)},localR{size_t(workers),std::vector<Agg>(N+1)};
  ranges(cap,workers,[&](size_t lo,size_t hi,int tid){auto&oneB=localB[tid];auto&oneR=localR[tid];for(uint32_t S=uint32_t(lo);S<uint32_t(hi);S++){int rank=__builtin_popcount(S);if(rank<3||!h[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t wb=q&-q,T=S^wb,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;int w=__builtin_ctz(wb);if(!(adj[w]&A)||!(adj[w]&B))continue;if(__builtin_ctz(A)>__builtin_ctz(B))std::swap(A,B);if(!(B&(uint32_t(1)<<r))||!(adj[w]&(uint32_t(1)<<r)))continue;int aa=__builtin_popcount(A),bb=__builtin_popcount(B);U term=h[S]*U(choose64(aa+bb,aa))*a[A]*ar[B];oneB[bb].img+=term;oneB[bb].states++;oneR[rank].img+=term;oneR[rank].states++;}}});
  for(int t=0;t<workers;t++)for(int k=1;k<=N;k++){byB[k].img+=localB[t][k].img;byRank[k].img+=localR[t][k].img;}
 }
 U image=0;for(auto z:byB)image+=z.img;
 static const U expectedA1[] = {0,0,4,29568ULL,35704394880ULL,U(8836673688265469760ULL)*10};
 U expectedA2 = n==5 ? U(7063261481362958280ULL)*1000 : (n==4 ? U(779299174080ULL) : U(143808ULL));
 if(a[full]!=expectedA1[n]||all.w!=expectedA2){std::cerr<<"global audit mismatch\n";return 4;}
 std::ofstream o(argv[2]);o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers<<",\n  \"state_array_bytes\": "<<uint64_t(3)*cap*sizeof(U)<<",\n  \"A1\": \""<<dec(a[full])<<"\",\n  \"A2_from_full_first_merge_scan\": \""<<dec(all.w)<<"\",\n  \"first_merge_state_records\": "<<all.states<<",\n  \"canonical_root_image_sum\": \""<<dec(image)<<"\",\n  \"by_canonical_B_size\": [";bool first=true;for(int k=1;k<=N;k++)if(byB[k].w){if(!first)o<<",";first=false;o<<"\n    {\"B_size\": "<<k<<", \"W\": \""<<dec(byB[k].w)<<"\", \"image\": \""<<dec(byB[k].img)<<"\", \"state_records\": "<<byB[k].states<<"}";}o<<"\n  ],\n  \"by_first_merge_rank\": [";first=true;for(int k=3;k<=N;k++)if(byRank[k].w){if(!first)o<<",";first=false;o<<"\n    {\"rank\": "<<k<<", \"W\": \""<<dec(byRank[k].w)<<"\", \"image\": \""<<dec(byRank[k].img)<<"\", \"state_records\": "<<byRank[k].states<<"}";}o<<"\n  ],\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";o.close();
 std::cout<<"n="<<n<<" A1="<<dec(a[full])<<" A2="<<dec(all.w)<<" image="<<dec(image)<<"\n";std::free(a);std::free(h);std::free(ar);
}
