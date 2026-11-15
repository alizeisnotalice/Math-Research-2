// Count beta>0: the lower peak's branch is a singleton at first merge.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
int main(int argc,char**argv){if(argc!=2)return 2;int n=std::stoi(argv[1]);if(n<2||n>5)return 2;int N=n*n;uint32_t size=1u<<N,full=size-1;std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(abs(v%n-u%n)+abs(v/n-u/n)==1)adj[v]|=1u<<u;
 U*a=(U*)calloc(size,sizeof(U)),*h=(U*)calloc(size,sizeof(U));if(!a||!h)return 3;
 for(uint32_t S=1;S<size;S++){if(!(S&(S-1))){a[S]=1;continue;}for(uint32_t bits=S;bits;bits&=bits-1){int v=__builtin_ctz(bits);uint32_t T=S^(1u<<v);if(adj[v]&T)a[S]+=a[T];}}
 h[full]=1;for(uint32_t S=full-1;S;S--)for(uint32_t bits=full^S;bits;bits&=bits-1){int v=__builtin_ctz(bits);if(adj[v]&S)h[S]+=h[S|(1u<<v)];}
 U support=0,both=0;for(uint32_t A=1;A<full;A++)if(a[A]){
  int k=__builtin_popcount(A);uint32_t closed=A;
  for(uint32_t bits=A;bits;bits&=bits-1)closed|=adj[__builtin_ctz(bits)];
  for(uint32_t bits=closed^A;bits;bits&=bits-1){int w=__builtin_ctz(bits);
   for(uint32_t cand=adj[w]&~closed;cand;cand&=cand-1){int v=__builtin_ctz(cand);U term=a[A]*k*h[A|(1u<<w)|(1u<<v)];support+=term;if(k==1)both+=term;}
  }
 }
 if(both*2!=a[full])return 4;
 std::cout<<"{\"n\":"<<n<<",\"A1\":"<<dec(a[full])<<",\"beta_positive_count\":"<<dec(support)<<",\"both_singleton_count\":"<<dec(both)<<",\"allocated_bytes\":"<<uint64_t(size)*32<<"}\n";free(a);free(h);
}
