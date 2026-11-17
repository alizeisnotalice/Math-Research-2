// Exact intersections of safe-root curved-cut events. Descending label order.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
static std::string dec(U x){if(!x)return "0";std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}std::reverse(s.begin(),s.end());return s;}
struct Constraint{uint32_t region,allowed,neighbors;};
int main(int argc,char**argv){
 if(argc<4){std::cerr<<"usage: N OUTPUT CURVE... (each curve is n digits, row left-widths)\n";return 2;}
 const int n=std::stoi(argv[1]),N=n*n;if(n<2||n>5)return 2;
 if(std::filesystem::exists(argv[2]))return 4;
 const uint32_t cap=uint32_t(1)<<N,full=cap-1;
 std::vector<uint32_t> adj(N);
 for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(std::abs(v%n-u%n)+std::abs(v/n-u/n)==1)adj[v]|=uint32_t(1)<<u;
 std::vector<std::vector<Constraint>> cs(N);
 std::vector<std::string> curves;
 for(int k=3;k<argc;k++){
  const std::string a=argv[k];if(int(a.size())!=n)return 2;
  uint32_t region[2]={0,0},allowed[2]={0,0};
  for(int y=0;y<n;y++)if(a[y]<'1'||a[y]>=char('0'+n))return 2;
  for(int v=0;v<N;v++)region[v%n>=a[v/n]-'0']|=uint32_t(1)<<v;
  for(int side=0;side<2;side++)for(int v=0;v<N;v++){
   const uint32_t bit=uint32_t(1)<<v;
   // Outer quarter bands, giving distance >= n-2*ceil(n/4)+1.
   const int band=(n+3)/4;
   const bool outer=side==0?v%n<band:v%n>=n-band;
   if((region[side]&bit)&&outer&&!(adj[v]&region[1-side]))allowed[side]|=bit;
  }
  for(int v=0;v<N;v++){const int side=(region[1]>>v)&1;cs[v].push_back({region[side],allowed[side],adj[v]&region[side]});}
  curves.push_back(a);
 }
 // Every count is <=25! < 2^84; no floating-point arithmetic.
 std::vector<U> dp(cap);dp[0]=1;uint64_t reached=0,transitions=0;
 for(uint32_t S=0;S<full;S++)if(dp[S]){
  ++reached;
  for(uint32_t bits=full^S;bits;bits&=bits-1){
   const uint32_t bit=bits&-bits;const int v=__builtin_ctz(bit);bool ok=true;
   for(const auto& c:cs[v])if((c.region&S)?!(c.neighbors&S):!(c.allowed&bit)){ok=false;break;}
   if(ok){dp[S|bit]+=dp[S];++transitions;}
  }
 }
 std::ofstream out(argv[2]);if(!out)return 5;
 out<<"{\"schema\":\"safe-curved-cut-intersection-v1\",\"n\":"<<n<<",\"curves\":[";
 for(size_t i=0;i<curves.size();i++){if(i)out<<',';out<<'\"'<<curves[i]<<'\"';}
 out<<"],\"count\":\""<<dec(dp[full])<<"\",\"reached_states\":"<<reached<<",\"valid_transitions\":"<<transitions<<",\"core_bytes\":"<<uint64_t(cap)*sizeof(U)<<"}\n";
 return 0;
}
