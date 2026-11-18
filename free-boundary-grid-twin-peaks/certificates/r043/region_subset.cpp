// Independent last-vertex subset DP, n7 profile domain with at most25 vertices.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using U=unsigned __int128;
static std::string dec(U x){if(!x)return "0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){
 if(argc!=3)return 2;const std::string h=argv[1];constexpr int n=7,N=49;
 if(h.size()!=7||std::filesystem::exists(argv[2]))return 2;
 std::vector<int> vs;std::array<int,N> idx;idx.fill(-1);
 for(int y=0;y<n;y++){if(h[y]!='3'&&h[y]!='4')return 2;for(int x=0;x<h[y]-'0';x++){idx[y*n+x]=int(vs.size());vs.push_back(y*n+x);}}
 const int m=int(vs.size());if(m>25){std::cerr<<"hard cap25 vertices (512MiB)\n";return 5;}
 std::vector<uint32_t> adj(m);std::vector<bool> allowed(m);
 for(int i=0;i<m;i++){
  const int v=vs[i],x=v%n,y=v/n;bool safe=true;
  for(auto [dx,dy]:std::array<std::pair<int,int>,4>{{{1,0},{-1,0},{0,1},{0,-1}}}){
   const int X=x+dx,Y=y+dy;if(X<0||X>=n||Y<0||Y>=n)continue;
   const int j=idx[Y*n+X];if(j<0)safe=false;else adj[i]|=uint32_t(1)<<j;
  }
  allowed[i]=x<2&&safe;
 }
 const uint32_t cap=uint32_t(1)<<m;std::vector<U> dp(cap);
 for(int i=0;i<m;i++)dp[uint32_t(1)<<i]=allowed[i]?1:0;
 for(uint32_t S=1;S<cap;S++){
  if(!(S&(S-1)))continue;U z=0;
  for(uint32_t q=S;q;q&=q-1){const uint32_t b=q&-q,T=S^b;const int i=__builtin_ctz(b);if(adj[i]&T)z+=dp[T];}
  dp[S]=z;
 }
 std::ofstream out(argv[2]);if(!out)return 4;
 out<<"{\"schema\":\"n7-region-last-vertex-cross-v1\",\"n\":7,\"profile\":\""<<h<<"\",\"vertices\":"<<m<<",\"rooted_count\":\""<<dec(dp.back())<<"\",\"core_bytes\":"<<uint64_t(cap)*sizeof(U)<<"}\n";
}
