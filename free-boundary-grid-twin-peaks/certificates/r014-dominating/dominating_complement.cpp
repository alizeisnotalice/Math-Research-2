// Exact coverage of the dominating-complement late-merge condition, n<=5.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
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
template<class F>static void parallelRanges(size_t c,int workers,F f){int used=std::min<int>(workers,int(c));std::vector<std::thread>ts;ts.reserve(used);for(int t=0;t<used;t++){size_t lo=c*size_t(t)/used,hi=c*size_t(t+1)/used;ts.emplace_back([=,&f](){f(lo,hi,t);});}for(auto&t:ts)t.join();}
static uint32_t component(uint32_t seed,uint32_t allowed,const std::array<uint32_t,25>&adj){uint32_t seen=seed,fr=seed;while(fr){uint32_t bit=fr&-fr;fr-=bit;uint32_t add=adj[__builtin_ctz(bit)]&allowed&~seen;seen|=add;fr|=add;}return seen;}
struct Agg{U w=0,d=0;uint64_t events=0,hfactFail=0;};
struct Part{U words=0;uint64_t states=0,hfactFail=0;};
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,workers=std::stoi(argv[3]);if(n<2||n>5||workers<1)return 2;
 uint32_t cap=uint32_t(1)<<N,full=cap-1;std::array<uint32_t,25>adj{};for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto alloc=[&](){return static_cast<U*>(std::calloc(cap,sizeof(U)));};U*a=alloc(),*mx=alloc(),*my=alloc(),*m2=alloc(),*h=alloc();if(!a||!mx||!my||!m2||!h){std::cerr<<"allocation failure\n";return 3;}auto start=std::chrono::steady_clock::now();
 for(int v=0;v<N;v++){uint32_t S=uint32_t(1)<<v;int x=v%n,y=v/n;a[S]=1;mx[S]=x;my[S]=y;m2[S]=U(x*x+y*y);}
 for(int k=2;k<=N;k++){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];U za=0,zx=0,zy=0,z2=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int v=__builtin_ctz(bit);uint32_t T=S^bit;if(adj[v]&T){za+=a[T];zx+=mx[T];zy+=my[T];z2+=m2[T];}}a[S]=za;mx[S]=zx;my[S]=zy;m2[S]=z2;}});}
 h[full]=1;for(int k=N-1;k>=1;k--){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i],out=full^S;U z=0;for(uint32_t q=out;q;q&=q-1){uint32_t bit=q&-q;if(adj[__builtin_ctz(bit)]&S)z+=h[S|bit];}h[S]=z;}});}
 std::array<U,26>fact{};fact[0]=1;for(int i=1;i<=N;i++)fact[i]=fact[i-1]*U(i);
 // category 0=only A, 1=only B, 2=both; partition key is (special-a, other-b, c0).
 int P=(N+1)*(N+1)*(N+1);auto pkey=[&](int aa,int bb,int c0){return (aa*(N+1)+bb)*(N+1)+c0;};
 std::vector<std::array<Agg,3>> local{size_t(workers)};std::vector<std::vector<Part>> plocal{size_t(workers),std::vector<Part>(P)};std::vector<std::array<Agg,2>> rlocal{size_t(workers)}; // ratio<1, ratio>=1
 for(int s=3;s<=N;s++){auto masks=masksOfSize(n,s);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int tid){auto&cats=local[tid];auto&parts=plocal[tid];auto&rat=rlocal[tid];for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];if(!h[S]||!a[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;uint32_t T=S^bit,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;uint32_t C=full^(A|B);bool pa=true,pb=true;for(uint32_t z=C;z;z&=z-1){uint32_t vb=z&-z;int v=__builtin_ctz(vb);if(!(adj[v]&A))pa=false;if(!(adj[v]&B))pb=false;}if(!pa&&!pb)continue;
    int as=__builtin_popcount(A),bs=__builtin_popcount(B),g=__builtin_popcount(C);U base=h[S]*U(choose64(s-1,as));U w=base*a[A]*a[B];U d=base*(a[B]*m2[A]+a[A]*m2[B]-U(2)*(mx[A]*mx[B]+my[A]*my[B]));int cat=pa&&pb?2:(pa?0:1);cats[cat].w+=w;cats[cat].d+=d;cats[cat].events++;
    // Choose the satisfying side; for both retain the canonical A (low-bit component).
    uint32_t X=A,Y=B;int xa=as,yb=bs;if(!pa&&pb){X=B;Y=A;xa=bs;yb=as;}
    int c0=0;for(uint32_t z=C;z;z&=z-1)if(adj[__builtin_ctz(z&-z)]&Y)c0++;if(c0<=0){std::cerr<<"unexpected c0=0\n";std::abort();}
    Part&pp=parts[pkey(xa,yb,c0)];pp.states++;pp.words+=w;if(h[S]!=fact[g-1]){pp.hfactFail++;cats[cat].hfactFail++;}
    // ratio=(c0/g)*C(as+bs,as)/C(N-1,yb-1); compare before division.
    uint64_t rn=uint64_t(c0)*choose64(as+bs,as),rd=uint64_t(g)*choose64(N-1,yb-1);Agg&rr=rat[rn<rd?0:1];rr.w+=w;rr.d+=d;rr.events++;
 }}});}
 std::array<Agg,3>cats{};std::array<Agg,2>rat{};std::vector<Part>parts(P);for(int t=0;t<workers;t++){for(int j=0;j<3;j++){cats[j].w+=local[t][j].w;cats[j].d+=local[t][j].d;cats[j].events+=local[t][j].events;cats[j].hfactFail+=local[t][j].hfactFail;}for(int j=0;j<2;j++){rat[j].w+=rlocal[t][j].w;rat[j].d+=rlocal[t][j].d;rat[j].events+=rlocal[t][j].events;}for(int k=0;k<P;k++){parts[k].words+=plocal[t][k].words;parts[k].states+=plocal[t][k].states;parts[k].hfactFail+=plocal[t][k].hfactFail;}}
 Agg cover;for(auto x:cats){cover.w+=x.w;cover.d+=x.d;cover.events+=x.events;cover.hfactFail+=x.hfactFail;}
 std::ofstream o(argv[2]);o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers<<",\n  \"state_array_bytes\": "<<uint64_t(5)*cap*sizeof(U)<<",\n  \"categories\": {";const char*names[3]={"only_A","only_B","both"};for(int j=0;j<3;j++){if(j)o<<", ";o<<"\""<<names[j]<<"\": {\"W\": \""<<dec(cats[j].w)<<"\", \"D\": \""<<dec(cats[j].d)<<"\", \"events\": "<<cats[j].events<<", \"h_factorial_failures\": "<<cats[j].hfactFail<<"}";}o<<"},\n  \"covered_union\": {\"W\": \""<<dec(cover.w)<<"\", \"D\": \""<<dec(cover.d)<<"\", \"events\": "<<cover.events<<", \"h_factorial_failures\": "<<cover.hfactFail<<"},\n  \"selected_side_ratio\": {\"less_than_one\": {\"W\": \""<<dec(rat[0].w)<<"\", \"D\": \""<<dec(rat[0].d)<<"\", \"events\": "<<rat[0].events<<"}, \"at_least_one\": {\"W\": \""<<dec(rat[1].w)<<"\", \"D\": \""<<dec(rat[1].d)<<"\", \"events\": "<<rat[1].events<<"}},\n  \"partition_layers\": [\n";bool first=true;for(int aa=1;aa<=N;aa++)for(int bb=1;bb<=N-aa;bb++)for(int c0=1;c0<=N-aa-bb;c0++){Part&p=parts[pkey(aa,bb,c0)];if(!p.states)continue;int g=N-aa-bb;if(p.states%uint64_t(c0)){std::cerr<<"state/c0 nondivisible\n";return 5;}uint64_t rn=uint64_t(c0)*choose64(aa+bb,aa),rd=uint64_t(g)*choose64(N-1,bb-1);if(!first)o<<",\n";first=false;o<<"    {\"special_a\": "<<aa<<", \"other_b\": "<<bb<<", \"g\": "<<g<<", \"c0\": "<<c0<<", \"state_count\": "<<p.states<<", \"partition_count\": "<<(p.states/c0)<<", \"family_word_count\": \""<<dec(p.words)<<"\", \"h_factorial_failures\": "<<p.hfactFail<<", \"lemma_ratio_numerator\": "<<rn<<", \"lemma_ratio_denominator\": "<<rd<<"}";}o<<"\n  ],\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";o.close();std::cout<<"n="<<n<<" coveredW="<<dec(cover.w)<<" coveredD="<<dec(cover.d)<<" events="<<cover.events<<"\n";std::free(a);std::free(mx);std::free(my);std::free(m2);std::free(h);
}
