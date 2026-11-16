// Exact first-merge mean/variance ledger with fixed-point rigorous bounds.
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
constexpr U L=U(1)<<18, L2=L*L;
static std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static uint64_t choose64(int n,int k){if(k>n-k)k=n-k;uint64_t z=1;for(int i=1;i<=k;i++)z=z*(n-k+i)/i;return z;}
static uint32_t nextComb(uint32_t x){uint32_t c=x&-x,r=x+c;return (((r^x)>>2)/c)|r;}
static std::vector<uint32_t> masksOfSize(int n,int k){int N=n*n;std::vector<uint32_t>z;z.reserve(size_t(choose64(N,k)));if(k==N){z.push_back((uint32_t(1)<<N)-1);return z;}uint32_t x=(uint32_t(1)<<k)-1,cap=uint32_t(1)<<N;while(x<cap){z.push_back(x);x=nextComb(x);}return z;}
template<class F>static void parallelRanges(size_t c,int workers,F f){int used=std::min<int>(workers,int(c));std::vector<std::thread>ts;ts.reserve(used);for(int t=0;t<used;t++){size_t lo=c*size_t(t)/used,hi=c*size_t(t+1)/used;ts.emplace_back([=,&f](){f(lo,hi,t);});}for(auto&t:ts)t.join();}
static uint32_t component(uint32_t seed,uint32_t allowed,const std::array<uint32_t,25>&adj){uint32_t seen=seed,fr=seed;while(fr){uint32_t bit=fr&-fr;fr-=bit;uint32_t add=adj[__builtin_ctz(bit)]&allowed&~seen;seen|=add;fr|=add;}return seen;}
struct Agg{U w=0,d=0,r2=0,meanLo=0,meanHi=0,withinLo=0,withinHi=0,betweenLo=0,betweenHi=0;};
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,workers=std::stoi(argv[3]);if(n<2||n>5||workers<1)return 2;
 U fact=1;for(int i=2;i<=N;i++)fact*=U(i);const U max=~U(0), bound=U(64)*L2;if(fact>max/bound){std::cerr<<"U128 width certificate fails\n";return 3;}
 uint32_t cap=uint32_t(1)<<N,full=cap-1;std::array<uint32_t,25>adj{};for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto alloc=[&](){return static_cast<U*>(std::calloc(cap,sizeof(U)));};U*a=alloc(),*mx=alloc(),*my=alloc(),*m2=alloc(),*h=alloc();if(!a||!mx||!my||!m2||!h){std::cerr<<"allocation failure\n";return 4;}auto start=std::chrono::steady_clock::now();
 for(int v=0;v<N;v++){uint32_t S=uint32_t(1)<<v;int x=v%n,y=v/n;a[S]=1;mx[S]=x;my[S]=y;m2[S]=U(x*x+y*y);}
 for(int k=2;k<=N;k++){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];U za=0,zx=0,zy=0,z2=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int v=__builtin_ctz(bit);uint32_t T=S^bit;if(adj[v]&T){za+=a[T];zx+=mx[T];zy+=my[T];z2+=m2[T];}}a[S]=za;mx[S]=zx;my[S]=zy;m2[S]=z2;}});}
 h[full]=1;for(int k=N-1;k>=1;k--){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i],out=full^S;U z=0;for(uint32_t q=out;q;q&=q-1){uint32_t bit=q&-q;if(adj[__builtin_ctz(bit)]&S)z+=h[S|bit];}h[S]=z;}});}
 const int MAXMIN=N/2;const size_t cells=size_t(N+1)*(MAXMIN+1);auto cell=[&](int s,int m){return size_t(s)*(MAXMIN+1)+m;};std::vector<std::vector<Agg>> local{size_t(workers),std::vector<Agg>(cells)};
 for(int s=3;s<=N;s++){auto masks=masksOfSize(n,s);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int tid){auto&ledger=local[tid];for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];if(!h[S]||!a[S])continue;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;uint32_t T=S^bit,A=component(T&-T,T,adj),B=T^A;if(!B||component(B&-B,B,adj)!=B)continue;int as=__builtin_popcount(A),bs=__builtin_popcount(B),mn=std::min(as,bs);U comb=choose64(s-1,as),base=h[S]*comb,aa=a[A],bb=a[B],w=base*aa*bb;
   U r2=base*(bb*m2[A]+aa*m2[B]);U dot=mx[A]*mx[B]+my[A]*my[B];U d=base*(bb*m2[A]+aa*m2[B]-U(2)*dot);
   auto fl=[&](U m,U den){return (L*m)/den;};auto ce=[&](U m,U den){return (L*m+den-1)/den;};
   U lxA=fl(mx[A],aa),lyA=fl(my[A],aa),lxB=fl(mx[B],bb),lyB=fl(my[B],bb);U uxA=ce(mx[A],aa),uyA=ce(my[A],aa),uxB=ce(mx[B],bb),uyB=ce(my[B],bb);
   U meanLo=w*(lxA*lxA+lyA*lyA+lxB*lxB+lyB*lyB);U meanHi=w*(uxA*uxA+uyA*uyA+uxB*uxB+uyB*uyB);
   U r2scaled=r2*L2,withinLo=(r2scaled>meanHi?r2scaled-meanHi:0),withinHi=r2scaled-meanLo;U dscaled=d*L2,betweenLo=(dscaled>withinHi?dscaled-withinHi:0),betweenHi=dscaled-withinLo;
   Agg&g=ledger[cell(s,mn)];g.w+=w;g.d+=d;g.r2+=r2;g.meanLo+=meanLo;g.meanHi+=meanHi;g.withinLo+=withinLo;g.withinHi+=withinHi;g.betweenLo+=betweenLo;g.betweenHi+=betweenHi;
 }}});}
 std::vector<Agg>ledger(cells);for(int t=0;t<workers;t++)for(size_t i=0;i<cells;i++){Agg&x=ledger[i];const Agg&y=local[t][i];x.w+=y.w;x.d+=y.d;x.r2+=y.r2;x.meanLo+=y.meanLo;x.meanHi+=y.meanHi;x.withinLo+=y.withinLo;x.withinHi+=y.withinHi;x.betweenLo+=y.betweenLo;x.betweenHi+=y.betweenHi;}
 Agg total;for(const auto&y:ledger){total.w+=y.w;total.d+=y.d;total.r2+=y.r2;total.meanLo+=y.meanLo;total.meanHi+=y.meanHi;total.withinLo+=y.withinLo;total.withinHi+=y.withinHi;total.betweenLo+=y.betweenLo;total.betweenHi+=y.betweenHi;}
 std::ofstream o(argv[2]);o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers<<",\n  \"fixed_point_L\": "<<dec(L)<<",\n  \"state_array_bytes\": "<<uint64_t(5)*cap*sizeof(U)<<",\n  \"u128_width_certificate\": \"64*L^2*"<<dec(fact)<<" < 2^128\",\n  \"totals\": {\"W\": \""<<dec(total.w)<<"\", \"D\": \""<<dec(total.d)<<"\", \"root_second_moment\": \""<<dec(total.r2)<<"\", \"mean_square_lower_numerator\": \""<<dec(total.meanLo)<<"\", \"mean_square_upper_numerator\": \""<<dec(total.meanHi)<<"\", \"within_lower_numerator\": \""<<dec(total.withinLo)<<"\", \"within_upper_numerator\": \""<<dec(total.withinHi)<<"\", \"between_lower_numerator\": \""<<dec(total.betweenLo)<<"\", \"between_upper_numerator\": \""<<dec(total.betweenHi)<<"\"},\n  \"rank_minbranch\": [\n";bool first=true;for(int s=3;s<=N;s++)for(int m=1;m<=MAXMIN;m++){const Agg&g=ledger[cell(s,m)];if(!g.w)continue;if(!first)o<<",\n";first=false;o<<"    {\"s\": "<<s<<", \"min_branch\": "<<m<<", \"W\": \""<<dec(g.w)<<"\", \"D\": \""<<dec(g.d)<<"\", \"root_second_moment\": \""<<dec(g.r2)<<"\", \"mean_square_lower_numerator\": \""<<dec(g.meanLo)<<"\", \"mean_square_upper_numerator\": \""<<dec(g.meanHi)<<"\", \"within_lower_numerator\": \""<<dec(g.withinLo)<<"\", \"within_upper_numerator\": \""<<dec(g.withinHi)<<"\", \"between_lower_numerator\": \""<<dec(g.betweenLo)<<"\", \"between_upper_numerator\": \""<<dec(g.betweenHi)<<"\"}";}o<<"\n  ],\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";o.close();std::cout<<"n="<<n<<" W="<<dec(total.w)<<" D="<<dec(total.d)<<"\n";std::free(a);std::free(mx);std::free(my);std::free(m2);std::free(h);
}
