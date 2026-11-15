// Count intersections of outer-half-root rectangular cut events by exact DP.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string dec(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return std::string(s.rbegin(),s.rend());}
struct Constraint {uint32_t region,allowed,neighbors;};
int main(int argc,char**argv){
 if(argc<3)return 2;int n=std::stoi(argv[1]),N=n*n;if(n<4||n>5)return 2;
 uint32_t size=1u<<N,full=size-1;std::vector<uint32_t>adj(N);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(abs(v%n-u%n)+abs(v/n-u/n)==1)adj[v]|=1u<<u;
 std::vector<std::vector<Constraint>>cs(N);
 for(int k=2;k<argc;k++){
  std::string a=argv[k];int axis=a[0]=='v'?0:1,j=std::stoi(a.substr(1));if(j<2||j>n-2)return 2;
  uint32_t regions[2]={0,0},allowed[2]={0,0};
  for(int v=0;v<N;v++){int x=axis==0?v%n:v/n;regions[x>=j]|=1u<<v;
    if(x<=(j-1)/2)allowed[0]|=1u<<v;if(x>=j+(n-j)/2)allowed[1]|=1u<<v;}
  for(int v=0;v<N;v++)for(int s=0;s<2;s++)if(regions[s]&(1u<<v))cs[v].push_back({regions[s],allowed[s],adj[v]&regions[s]});
 }
 U*f=(U*)calloc(size,sizeof(U));if(!f)return 3;f[0]=1;
 for(uint32_t S=0;S<full;S++)if(f[S])for(uint32_t bits=full^S;bits;bits&=bits-1){
  int v=__builtin_ctz(bits);uint32_t b=1u<<v;bool ok=true;
  for(const auto&c:cs[v])if((c.region&S)?!(c.neighbors&S):!(c.allowed&b)){ok=false;break;}
  if(ok)f[S|b]+=f[S];
 }
 std::cout<<dec(f[full])<<"\n";free(f);
}
