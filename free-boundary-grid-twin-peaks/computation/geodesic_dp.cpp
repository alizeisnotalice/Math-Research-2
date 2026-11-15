#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
int n,N;uint32_t full;std::vector<uint32_t>adj;std::vector<U>num;std::vector<uint64_t>paths;U*h;
void extend(int root,int v,uint32_t S,int d){num[d]+=h[S];paths[d]++;for(uint32_t bits=adj[v]&~S;bits;bits&=bits-1){int u=__builtin_ctz(bits);if(abs(u%n-root%n)+abs(u/n-root/n)==d)extend(root,u,S|(1u<<u),d+1);}}
int main(int argc,char**argv){if(argc!=2)return 2;n=std::stoi(argv[1]);if(n<2||n>5)return 2;N=n*n;full=(1u<<N)-1;adj.resize(N);num.resize(2*n);paths.resize(2*n);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(abs(v%n-u%n)+abs(v/n-u/n)==1)adj[v]|=1u<<u;
 h=(U*)calloc(size_t(full)+1,sizeof(U));if(!h)return 3;h[full]=1;
 for(uint32_t S=full-1;S;S--)for(uint32_t bits=full^S;bits;bits&=bits-1){int v=__builtin_ctz(bits);if(adj[v]&S)h[S]+=h[S|(1u<<v)];}
 for(int root=0;root<N;root++)extend(root,root,1u<<root,1);
 std::cout<<"{\"n\":"<<n<<",\"A1\":"<<dec(num[1])<<",\"rows\":[";
 for(int d=1;d<2*n;d++){if(d>1)std::cout<<",";std::cout<<"{\"d\":"<<d<<",\"numerator\":"<<dec(num[d])<<",\"ordered_prefixes\":"<<paths[d]<<"}";}
 std::cout<<"]}\n";free(h);
}
