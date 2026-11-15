// Exact count for a prescribed unordered peak pair, independent of merge decomposition.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
int main(int argc,char**argv){
 if(argc!=4)return 2;int n=std::stoi(argv[1]),N=n*n,p=std::stoi(argv[2]),q=std::stoi(argv[3]);
 if(n<2||n>5||p<0||q<0||p>=N||q>=N||p==q)return 2;
 uint32_t size=1u<<N,full=size-1,roots=(1u<<p)|(1u<<q);std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(abs(v%n-u%n)+abs(v/n-u/n)==1)adj[v]|=1u<<u;
 U*f=(U*)calloc(size,sizeof(U));if(!f)return 3;f[0]=1;
 for(uint32_t S=0;S<full;S++)if(f[S])for(uint32_t bits=full^S;bits;bits&=bits-1){
  int v=__builtin_ctz(bits);uint32_t b=1u<<v;
  if(bool(roots&b)!=bool(adj[v]&S))f[S|b]+=f[S];
 }
 std::cout<<dec(f[full])<<"\n";free(f);
}
