// Exact n<=4 audit of the early-interface completion-gain lemma.
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
 #include <gmpxx.h>
using U=__uint128_t; using Big=mpz_class;
static std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
static std::string decbig(const Big& x){return x.get_str();}
static Big big(U x){Big z=static_cast<unsigned long>(x>>64);z<<=64;z+=static_cast<unsigned long>(x);return z;}
static uint64_t choose64(int n,int k){if(k>n-k)k=n-k;uint64_t z=1;for(int i=1;i<=k;i++)z=z*(n-k+i)/i;return z;}
static uint32_t nextComb(uint32_t x){uint32_t c=x&-x,r=x+c;return (((r^x)>>2)/c)|r;}
static std::vector<uint32_t> masksOfSize(int n,int k){int N=n*n;std::vector<uint32_t>z;z.reserve(size_t(choose64(N,k)));if(k==N){z.push_back((uint32_t(1)<<N)-1);return z;}uint32_t x=(uint32_t(1)<<k)-1,cap=uint32_t(1)<<N;while(x<cap){z.push_back(x);x=nextComb(x);}return z;}
template<class F>static void parallelRanges(size_t c,int workers,F f){int used=std::min<int>(workers,int(c));std::vector<std::thread>ts;ts.reserve(used);for(int t=0;t<used;t++){size_t lo=c*size_t(t)/used,hi=c*size_t(t+1)/used;ts.emplace_back([=,&f](){f(lo,hi,t);});}for(auto&t:ts)t.join();}
static uint32_t component(uint32_t seed,uint32_t allowed,const std::array<uint32_t,25>&adj){uint32_t seen=seed,fr=seed;while(fr){uint32_t b=fr&-fr;fr-=b;uint32_t add=adj[__builtin_ctz(b)]&allowed&~seen;seen|=add;fr|=add;}return seen;}
struct Dir{uint32_t lo,hi;U hside;bool sideLo;}; struct Agg{U w=0,d=0;uint64_t pairs=0;};
struct GammaEx{bool set=false;U f=0,h=1;uint32_t A=0,B=0;bool sideA=true;};
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,workers=std::stoi(argv[3]);if(n<2||n>4||workers<1)return 2;
 uint32_t cap=uint32_t(1)<<N,full=cap-1;std::array<uint32_t,25>adj{};for(int y=0;y<n;y++)for(int x=0;x<n;x++){int v=y*n+x;for(auto d:std::array<std::array<int,2>,4>{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}}}}){int X=x+d[0],Y=y+d[1];if(X>=0&&X<n&&Y>=0&&Y<n)adj[v]|=uint32_t(1)<<(Y*n+X);}}
 auto alloc=[&](){return static_cast<U*>(std::calloc(cap,sizeof(U)));};U*a=alloc(),*mx=alloc(),*my=alloc(),*m2=alloc(),*h=alloc();if(!a||!mx||!my||!m2||!h)return 3;auto start=std::chrono::steady_clock::now();
 for(int v=0;v<N;v++){uint32_t S=uint32_t(1)<<v;int x=v%n,y=v/n;a[S]=1;mx[S]=x;my[S]=y;m2[S]=U(x*x+y*y);}
 for(int k=2;k<=N;k++){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i];U za=0,zx=0,zy=0,z2=0;for(uint32_t q=S;q;q&=q-1){uint32_t bit=q&-q;int v=__builtin_ctz(bit);uint32_t T=S^bit;if(adj[v]&T){za+=a[T];zx+=mx[T];zy+=my[T];z2+=m2[T];}}a[S]=za;mx[S]=zx;my[S]=zy;m2[S]=z2;}});}
 h[full]=1;for(int k=N-1;k>=1;k--){auto masks=masksOfSize(n,k);parallelRanges(masks.size(),workers,[&](size_t lo,size_t hi,int){for(size_t i=lo;i<hi;i++){uint32_t S=masks[i],out=full^S;U z=0;for(uint32_t q=out;q;q&=q-1){uint32_t bit=q&-q;if(adj[__builtin_ctz(bit)]&S)z+=h[S|bit];}h[S]=z;}});}
 std::vector<uint32_t> conn;for(uint32_t B=1;B<cap;B++)if(a[B])conn.push_back(B);
 std::vector<std::vector<Dir>> local{size_t(workers)};
 parallelRanges(conn.size(),workers,[&](size_t lo,size_t hi,int tid){U*hb=static_cast<U*>(std::calloc(cap,sizeof(U)));auto&out=local[tid];for(size_t j=lo;j<hi;j++){uint32_t B=conn[j],H=full^B,nb=0;for(uint32_t q=B;q;q&=q-1)nb|=adj[__builtin_ctz(q&-q)];uint32_t eligible=H&~nb;std::memset(hb,0,size_t(cap)*sizeof(U));hb[H]=1;for(uint32_t S=H;;S=(S-1)&H){if(S!=H){U z=0,outside=H^S;for(uint32_t q=outside;q;q&=q-1){uint32_t bit=q&-q;if(adj[__builtin_ctz(bit)]&S)z+=hb[S|bit];}hb[S]=z;}if(!S)break;}
   for(uint32_t A=eligible;A;A=(A-1)&eligible){if(!a[A])continue;uint32_t C=full^(A|B);U f=0;for(uint32_t q=C;q;q&=q-1){uint32_t bit=q&-q;int w=__builtin_ctz(bit);if((adj[w]&A)&&(adj[w]&B))f+=h[A|B|bit];}if(!f)continue;bool sideLo=__builtin_ctz(A)<__builtin_ctz(B);out.push_back({sideLo?A:B,sideLo?B:A,hb[A],sideLo});}
 }std::free(hb);});
 std::vector<Dir> dirs;for(auto&v:local)dirs.insert(dirs.end(),v.begin(),v.end());std::sort(dirs.begin(),dirs.end(),[](const Dir&x,const Dir&y){uint64_t kx=(uint64_t(x.lo)<<32)|x.hi,ky=(uint64_t(y.lo)<<32)|y.hi;if(kx!=ky)return kx<ky;return x.sideLo<y.sideLo;});
 std::array<Agg,3>avail{}; // only lower-A, only upper-B, both
 Agg allpairs{},unavailable{},pos1{},pos2{},pos4{},neg2{},neg4{};Big sumUpper=0;const Big scale=Big(1)<<128;std::vector<Agg> strata(size_t(N+1)*(N/2+1));auto ix=[&](int s,int m){return size_t(s)*(N/2+1)+m;};uint64_t records=0;GammaEx gmin,gmax,zeroExample;
 for(size_t i=0;i<dirs.size();){size_t j=i+1;while(j<dirs.size()&&dirs[j].lo==dirs[i].lo&&dirs[j].hi==dirs[i].hi)j++;if(j-i!=2){std::cerr<<"missing reverse direction\n";return 4;}U h0=dirs[i].sideLo?dirs[i].hside:dirs[i+1].hside,h1=dirs[i].sideLo?dirs[i+1].hside:dirs[i].hside;uint32_t A=dirs[i].lo,B=dirs[i].hi,C=full^(A|B);int aa=__builtin_popcount(A),bb=__builtin_popcount(B),g=__builtin_popcount(C),s=aa+bb+1,mn=std::min(aa,bb);U f=0;for(uint32_t q=C;q;q&=q-1){uint32_t bit=q&-q;int w=__builtin_ctz(bit);if((adj[w]&A)&&(adj[w]&B))f+=h[A|B|bit];}U base=U(choose64(aa+bb,aa))*f,w=base*a[A]*a[B],d=base*(a[B]*m2[A]+a[A]*m2[B]-U(2)*(mx[A]*mx[B]+my[A]*my[B]));allpairs.w+=w;allpairs.d+=d;allpairs.pairs++;bool u0=h0>0,u1=h1>0;if(!u0&&!u1){unavailable.w+=w;unavailable.d+=d;unavailable.pairs++;if(!zeroExample.set)zeroExample={true,f,0,A,B,true};i=j;continue;}int cat=u0&&u1?2:(u0?0:1);avail[cat].w+=w;avail[cat].d+=d;avail[cat].pairs++;strata[ix(s,mn)].w+=w;strata[ix(s,mn)].d+=d;strata[ix(s,mn)].pairs++;records++;
  auto note=[&](U hh,bool sideA){GammaEx z{true,f,hh,A,B,sideA};if(!gmin.set||big(f)*big(gmin.h)<big(gmin.f)*big(hh))gmin=z;if(!gmax.set||big(f)*big(gmax.h)>big(gmax.f)*big(hh))gmax=z;};if(u0)note(h0,true);if(u1)note(h1,false);
  bool p1=(u0&&f<=h0)||(u1&&f<=h1);bool p2=(u0&&f<=(h0<<g))||(u1&&f<=(h1<<g));bool p4=(u0&&f<=(h0<<(2*g)))||(u1&&f<=(h1<<(2*g)));bool t2=(u0&&(f<<g)<=h0)||(u1&&(f<<g)<=h1);bool t4=(u0&&(f<<(2*g))<=h0)||(u1&&(f<<(2*g))<=h1);if(p1){pos1.w+=w;pos1.d+=d;pos1.pairs++;}if(p2){pos2.w+=w;pos2.d+=d;pos2.pairs++;}if(p4){pos4.w+=w;pos4.d+=d;pos4.pairs++;}if(t2){neg2.w+=w;neg2.d+=d;neg2.pairs++;}if(t4){neg4.w+=w;neg4.d+=d;neg4.pairs++;}
  U q0=0,q1=0;if(u0)q0=a[A]*a[B]*h0*U(choose64(N,bb)-choose64(bb+g-1,bb));if(u1)q1=a[A]*a[B]*h1*U(choose64(N,aa)-choose64(aa+g-1,aa));U dstar=std::max(q0,q1);Big numerator=big(w)<<128,denom=big(dstar);sumUpper+=(numerator+denom-1)/denom;i=j;
 }
 Agg uni;for(auto z:avail){uni.w+=z.w;uni.d+=z.d;uni.pairs+=z.pairs;}
 Big sumLower=sumUpper-Big(static_cast<unsigned long>(records));
 auto emitAgg=[](std::ostream&z,const Agg&x){z<<"{\"W\": \""<<dec(x.w)<<"\", \"D\": \""<<dec(x.d)<<"\", \"pairs\": "<<x.pairs<<"}";};
 auto emitEx=[](std::ostream&z,const GammaEx&x){z<<"{\"set\": "<<(x.set?"true":"false");if(x.set)z<<", \"A_mask\": "<<x.A<<", \"B_mask\": "<<x.B<<", \"special_side\": \""<<(x.sideA?"canonical_A":"canonical_B")<<"\", \"Fnum\": \""<<dec(x.f)<<"\", \"H\": \""<<dec(x.h)<<"\"";z<<"}";};
 std::ofstream o(argv[2]);o<<"{\n  \"n\": "<<n<<",\n  \"workers\": "<<workers<<",\n  \"global_state_array_bytes\": "<<uint64_t(5)*cap*sizeof(U)<<",\n  \"connected_B\": "<<conn.size()<<",\n  \"ordered_direction_records\": "<<dirs.size()<<",\n  \"all_pairs\": ";emitAgg(o,allpairs);o<<",\n  \"unavailable_both_H_zero\": ";emitAgg(o,unavailable);o<<",\n  \"available_union\": ";emitAgg(o,uni);o<<",\n  \"availability\": {\"only_canonical_A\": ";emitAgg(o,avail[0]);o<<", \"only_canonical_B\": ";emitAgg(o,avail[1]);o<<", \"both\": ";emitAgg(o,avail[2]);o<<"},\n  \"gamma_union_positive_exponents\": {\"le_1\": ";emitAgg(o,pos1);o<<", \"le_2_pow_g\": ";emitAgg(o,pos2);o<<", \"le_4_pow_g\": ";emitAgg(o,pos4);o<<"},\n  \"gamma_union_negative_exponents\": {\"le_2_neg_g\": ";emitAgg(o,neg2);o<<", \"le_4_neg_g\": ";emitAgg(o,neg4);o<<"},\n  \"gamma_examples\": {\"minimum\": ";emitEx(o,gmin);o<<", \"maximum\": ";emitEx(o,gmax);o<<", \"both_H_zero\": ";emitEx(o,zeroExample);o<<"},\n  \"sum_F_over_Dstar_rigorous_interval\": {\"lower_numerator\": \""<<decbig(sumLower)<<"\", \"upper_numerator\": \""<<decbig(sumUpper)<<"\", \"denominator\": \"340282366920938463463374607431768211456\", \"rounding\": \"termwise ceiling on 2^128 grid; lower is upper minus usable-pair count\"},\n  \"rank_minbranch\": [\n";bool first=true;for(int ss=3;ss<=N;ss++)for(int mm=1;mm<=N/2;mm++){const Agg&z=strata[ix(ss,mm)];if(!z.w)continue;if(!first)o<<",\n";first=false;o<<"    {\"s\": "<<ss<<", \"min_branch\": "<<mm<<", \"W\": \""<<dec(z.w)<<"\", \"D\": \""<<dec(z.d)<<"\", \"pairs\": "<<z.pairs<<"}";}o<<"\n  ],\n  \"wall_seconds\": "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";o.close();std::cout<<"n="<<n<<" availableW="<<dec(uni.w)<<" D="<<dec(uni.d)<<" pairs="<<uni.pairs<<"\n";std::free(a);std::free(mx);std::free(my);std::free(m2);std::free(h);
}
