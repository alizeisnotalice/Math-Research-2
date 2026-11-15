// Exact h({u,v}) table for n<=5, where h(S) is the number of descending
// growth orders after every point of S has already been revealed.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using U = unsigned __int128;
static std::string out(U x){if(!x)return"0";std::string s;while(x){s.push_back('0'+x%10);x/=10;}return {s.rbegin(),s.rend()};}
int main(int argc,char**argv){
  if(argc!=2)return 2;int n=std::atoi(argv[1]);if(n<2||n>5)return 2;int N=n*n;uint32_t full=(1u<<N)-1;
  uint32_t col0=0,coll=0;for(int y=0;y<n;y++){col0|=1u<<(y*n);coll|=1u<<(y*n+n-1);}
  auto neighbours=[&](uint32_t s){return (((s&~col0)>>1)|((s&~coll)<<1)|(s<<n)|(s>>n))&full;};
  std::vector<U> h(uint64_t(full)+1);h[full]=1;
  for(uint32_t s=full-1;s;s--){uint32_t q=neighbours(s)&~s;U z=0;while(q){uint32_t b=q&-q;q-=b;z+=h[s|b];}h[s]=z;}
  int diam=2*(n-1);std::vector<U> by(diam+1);std::vector<uint64_t> pairs(diam+1);U total=0;
  for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){
    int d=std::abs(u%n-v%n)+std::abs(u/n-v/n);if(d==1)continue;U z=h[(1u<<u)|(1u<<v)];by[d]+=z;pairs[d]++;total+=z;
  }
  std::array<std::pair<const char*,uint32_t>,3> probes{{
    {"opposite_corners",(1u<<0)|(1u<<(N-1))},
    {"middle_straight_distance2",(1u<<(N/2))|(1u<<(N/2+2))},
    {"top_side_span",(1u<<0)|(1u<<(n-1))}
  }};
  std::cout<<"{\"n\":"<<n<<",\"unordered_top2_sum\":\""<<out(total)<<"\",\"ordered_top2_sum\":\""<<out(2*total)<<"\",\"by_distance\":{";
  bool first=true;for(int d=2;d<=diam;d++)if(pairs[d]){if(!first)std::cout<<',';first=false;std::cout<<'\"'<<d<<"\":{\"pairs\":"<<pairs[d]<<",\"h_sum\":\""<<out(by[d])<<"\"}";}
  std::cout<<"},\"selected_h\":{";for(size_t i=0;i<probes.size();i++){if(i)std::cout<<',';std::cout<<'\"'<<probes[i].first<<"\":\""<<out(h[probes[i].second])<<"\"";}std::cout<<"}}\n";
}
