// Exact total normalized promotion mass, via first touch in one-peak growth.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
int main(int argc,char**argv){if(argc!=2)return 2;int n=std::stoi(argv[1]);if(n<2||n>5)return 2;int N=n*n;uint32_t size=1u<<N,full=size-1;std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(abs(v%n-u%n)+abs(v/n-u/n)==1)adj[v]|=1u<<u;
 U*a=(U*)calloc(size,sizeof(U)),*h=(U*)calloc(size,sizeof(U));if(!a||!h)return 3;
 for(uint32_t S=1;S<size;S++){if(!(S&(S-1))){a[S]=1;continue;}for(uint32_t bits=S;bits;bits&=bits-1){int v=__builtin_ctz(bits);uint32_t T=S^(1u<<v);if(adj[v]&T)a[S]+=a[T];}}
 h[full]=1;for(uint32_t S=full-1;S;S--)for(uint32_t bits=full^S;bits;bits&=bits-1){int v=__builtin_ctz(bits);if(adj[v]&S)h[S]+=h[S|(1u<<v)];}
 uint64_t D=1;for(int i=1;i<=N-2;i++)D=std::lcm(D,uint64_t(i));std::vector<U>omega(N);
 for(int k=1;k<N-1;k++)for(int t=N-k;t<=N-1;t++)omega[k]+=D/(t-1);
 U mass=0;for(uint32_t S=1;S<full;S++)if(a[S]){
  int k=__builtin_popcount(S);if(k>=N-1)continue;uint32_t closed=S;
  for(uint32_t bits=S;bits;bits&=bits-1)closed|=adj[__builtin_ctz(bits)];
  uint32_t front=closed^S;
  for(uint32_t bits=front;bits;bits&=bits-1){int w=__builtin_ctz(bits);int newtouch=__builtin_popcount(adj[w]&~closed);
   if(newtouch)mass+=a[S]*h[S|(1u<<w)]*newtouch*omega[k];
  }
 }
 std::cout<<"{\"n\":"<<n<<",\"A1\":"<<dec(a[full])<<",\"mass_numerator\":"<<dec(mass)<<",\"mass_denominator\":"<<D<<",\"allocated_bytes\":"<<uint64_t(size)*32<<"}\n";free(a);free(h);
}
