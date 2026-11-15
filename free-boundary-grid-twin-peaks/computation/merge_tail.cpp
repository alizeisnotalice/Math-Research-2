// Exact threshold tail dynamic program for the first-merger rank.
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

struct __attribute__((packed)) U80 { uint64_t lo; uint16_t hi; };
static inline U80 add(U80 a,U80 b) { U80 c{uint64_t(a.lo+b.lo),uint16_t(a.hi+b.hi+(a.lo+b.lo<a.lo))}; return c; }
static inline bool zero(U80 a){return !(a.lo||a.hi);}
static inline U80 mulsmall(U80 a,uint32_t b) { unsigned __int128 z=(unsigned __int128)a.lo*b; return {uint64_t(z),uint16_t((z>>64)+a.hi*b)}; }
static inline U80 mul(U80 a,U80 b) { unsigned __int128 z=(unsigned __int128)a.lo*b.lo; return {uint64_t(z),uint16_t((z>>64)+a.lo*b.hi+a.hi*b.lo)}; }
static std::string show(U80 x){if(zero(x))return"0";std::string s;while(!zero(x)){uint64_t carry=x.hi%10; uint16_t nhi=x.hi/10; unsigned __int128 z=((unsigned __int128)carry<<64)|x.lo; uint64_t nlo=(uint64_t)(z/10), rem=(uint64_t)(z%10); x={nlo,nhi};s.push_back('0'+rem);}return std::string(s.rbegin(),s.rend());}
static uint32_t adj[25]; static int n,N; static uint32_t ALL;
static uint32_t comp(uint32_t S){uint32_t seen=0,todo=S&-S;while(todo){uint32_t b=todo&-todo;todo^=b;int v=__builtin_ctz(b);seen|=b;todo|=adj[v]&S&~seen;}return seen;}
static uint32_t choose[26][26];
static U80 disconnected_weight(uint32_t S,const U80*a){if(!S)return{};uint32_t A=comp(S),B=S^A;if(!B||comp(B)!=B)return{};int m=__builtin_popcount(S),r=__builtin_popcount(A);return mulsmall(mul(a[A],a[B]),choose[m][r]);}
int main(int argc,char**argv){if(argc!=3)return 2;n=std::stoi(argv[1]);int threshold=std::stoi(argv[2]);N=n*n;if(n<1||n>5)return 2;ALL=(1u<<N)-1;for(int i=0;i<=25;i++){choose[i][0]=choose[i][i]=1;for(int j=1;j<i;j++)choose[i][j]=choose[i-1][j-1]+choose[i-1][j];}
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=1u<<u;
 size_t count=(size_t)ALL+1;U80*a=(U80*)std::calloc(count,sizeof(U80));U80*g=(U80*)std::calloc(count,sizeof(U80));if(!a||!g)return 3;auto start=std::chrono::steady_clock::now();
 for(uint32_t S=1;S<=ALL;S++){int m=__builtin_popcount(S);if(!(S&(S-1))){a[S]={1,0};continue;}for(uint32_t z=S;z;z&=z-1){int v=__builtin_ctz(z);uint32_t T=S^(1u<<v);if(adj[v]&T)a[S]=add(a[S],a[T]);}
  // g counts two-birth orders whose first merger occurs strictly after threshold.
  if(comp(S)==S) for(uint32_t z=S;z;z&=z-1){int v=__builtin_ctz(z);uint32_t T=S^(1u<<v);if(!T)continue;uint32_t A=comp(T),B=T^A;if(!B){if(adj[v]&T)g[S]=add(g[S],g[T]);}else if(comp(B)==B && (adj[v]&A) && (adj[v]&B) && m>threshold)g[S]=add(g[S],disconnected_weight(T,a));}
 }
 std::cout<<"{\"n\":"<<n<<",\"threshold\":"<<threshold<<",\"tail_count\":\""<<show(g[ALL])<<"\",\"A1\":\""<<show(a[ALL])<<"\",\"allocated_bytes\":"<<2*count*sizeof(U80)<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}
